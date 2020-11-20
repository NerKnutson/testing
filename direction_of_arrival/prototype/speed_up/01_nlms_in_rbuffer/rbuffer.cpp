#include "rbuffer.h"

void RingBuffer::Reset()
{
	//std::lock_guard<std::mutex> lock(mutex_);
	head_ = tail_;
	full_ = false;
	//std::lock_guard<std::mutex> unlock(mutex_);
}

// State Tracking
bool RingBuffer::Empty() const
{
// if head and tail are equal, we are empty
	return (!full_ && (head_ == tail_));
}

bool RingBuffer::Full() const
{
// if tail is ahead the head by 1, we are full
	return full_;
}

size_t RingBuffer::Capacity() const
{
// returns the capacity of the ring buffer
	return max_size_;
}

size_t RingBuffer::Size() const
{
// returns the current size of the buffer
	size_t size = max_size_;
	if(!full_)
	{
		if(head_ >= tail_)
		{
			size = head_ - tail_;
		}
		else
		{
			size = max_size_ + head_ - tail_;
		}
	}
	return size;
}

// Adding Data
void RingBuffer::Put(const double item[])
{
// places the new data item[] at the head of the buffer and increments the head
	//std::lock_guard<std::mutex> lock(mutex_);
	for(int i = 0; i < m_buffer_number; i++)
	{
		buf_[i][head_] = item[i];
	}
	if(full_)
	{
		tail_ = (tail_ + 1) % max_size_;
	}
	head_ = (head_ + 1) % max_size_;
	full_ = head_ == tail_;
	//std::lock_guard<std::mutex> unlock(mutex_);
}
void RingBuffer::UnPut(int size)
{
// Soft-Removes by placing the head of the buffer back by size
	if(size >= Size())
		Reset();
	else
		head_ = (head_ - size) % max_size_;
}

// Retrieving Data
void RingBuffer::Get(double returned_data[])
{
// retrieves data from the tail and increments the tail (leaving that entry to be rewritten)
	//std::lock_guard<std::mutex> lock(mutex_);

	if(Empty())
	{
		returned_data = NULL;
	}
// Read data and advance the tail (we now have a free space)
	for(int i = 0; i < m_buffer_number; i++)
	{
		returned_data[i] =  buf_[i][tail_];
	}
	full_ = false;
	tail_ = (tail_ + 1) % max_size_;
	//std::lock_guard<std::mutex> unlock(mutex_);
}

// Looks at the index from the tail and returns data for a particular channel
// does not increment tail (essentially leaves buffer alone)
double RingBuffer::Look(int index, int channel)
{
	//std::lock_guard<std::mutex> lock(mutex_);
	if(Empty())
	{
		return (double) NULL;
	}
	// Return buffer for <channel> at tail + <index>
	//std::lock_guard<std::mutex> unlock(mutex_);
	return buf_[channel][(tail_+index) % max_size_];
}

// Using the Look function, dumps all data to history[N_chan][marker] and shot[N_chan][max_size_ - marker]
void RingBuffer::Dump(int marker, double* history, double* shot)
{
	//std::lock_guard<std::mutex> lock(mutex_);
	if(Empty())
	{
		history = NULL;
		shot = NULL;
	}
// Read data
	for(int i = 0; i < marker; i++)
	{
		for(int j = 0; j < m_buffer_number; j++)
			history[i + j*marker] = Look(i, j);
	}
	for(int i = marker; i < max_size_; i++)
	{
		for(int j = 0; j < m_buffer_number; j++)
			shot[i-marker + j*(max_size_-marker)] = Look(i, j);
	}
	//std::lock_guard<std::mutex> unlock(mutex_);
}


bool RingBuffer::Process(double data[])
{
	if (index < m_filter_order)
	{
	// Fill array with data until we have enough to begin NLMS
		for (int i = 0; i < m_buffer_number; i++)
		{
			m_x[i][index] = data[i];
		}
		index++;
		return false;
	}
	else
	{
		index++;
		index = m_filter_order + index%m_filter_order;
	// Process m_x with present data[]
		double place_holder = m_tiny;
// Continuous Learning
		for (int i = 0; i < m_buffer_number; i++)
		{
			m_error[i] = data[i];
			place_holder += data[i]*data[i];
		}
	// ERROR -= OLD_DATA
		for (int i = 0; i < m_filter_order; i++)
		{
			for (int j = 0; j < m_buffer_number; j++)
			{
				place_holder += m_x[j][i]*m_x[j][i];
				for (int k = 0; k < m_buffer_number; k++)
				{
					m_error[j] -= m_x[k][(i+index)%m_filter_order]*m_coeff[j][k][(i+index)%m_filter_order];
					if (j != k)
						m_error[j] -= data[k]*m_coeff[j][k][m_filter_order];
				}
			}
		}
	// Calculating new coefficients
	// new_coeff = old_coeff + error*data*step_size/|data|^2
		place_holder = m_step_size/place_holder;
		for (int i = 0; i < m_filter_order; i++)
		{
			for (int j = 0; j < m_buffer_number; j++)
				for (int k = 0; k < m_buffer_number; k++)
				{
					if(this->IsRecording())
					{
						 m_coeff[j][k][(i+index)%m_filter_order] += place_holder*m_error[j]*m_x[k][(i+index)%m_filter_order];
					 if (i == m_step_size - 1 and j != k)
						 m_coeff[j][k][m_filter_order] += place_holder*m_error[j]*data[k];
					}
				}
		}
	// Check for shot and increment variance
		for (int i = 0; i < m_buffer_number; i++)
		{
			m_record_history[i] = m_triggered[i] = (m_error[i]*m_error[i] > m_threshold*m_var[i]);
			if (!m_record_history[i] and !m_var_recorded[i])
			{
				m_triggered_var[i] = m_var[i];
				m_var_recorded[i] = true;
			}
			if (m_var[i] <= 1.20*m_triggered_var[i])
			{
				m_record_history[i] = true;
				m_var_recorded[i] = false;
			}

			m_var[i] = m_alpha*m_var[i] + (1-m_alpha)*m_error[i]*m_error[i];
		// Add data to m_x for next round; If triggered, does not add that data point
			m_x[i][index%m_filter_order] = data[i]*m_triggered[i];
		}
		return true;
	}
		//  ERROR = PRESENT_DATA
/*
		for (int i = 0; i < m_buffer_number; i++)
		{
			m_error[i] = data[i]*(!m_triggered[i]);
			place_holder += data[i]*data[i]*(!m_triggered[i]);
		}
	// ERROR -= OLD_DATA
		for (int i = 0; i < m_filter_order; i++)
		{
			for (int j = 0; j < m_buffer_number; j++)
			{
				place_holder += m_x[j][i]*m_x[j][i]*(!m_triggered[j]);
				for (int k = 0; k < m_buffer_number; k++)
				{
					m_error[j] -= m_x[k][(i+index)%m_filter_order]*m_coeff[j][k][(i+index)%m_filter_order]*(!m_triggered[j]);
					if (j != k)
						m_error[j] -= data[k]*m_coeff[j][k][m_filter_order]*(!m_triggered[j]);
				}
			}
		}
	// Calculating new coefficients
	// new_coeff = old_coeff + error*data*step_size/|data|^2
		place_holder = m_step_size/place_holder;
		for (int i = 0; i < m_filter_order; i++)
		{
			for (int j = 0; j < m_buffer_number; j++)
				for (int k = 0; k < m_buffer_number; k++)
				{
					 m_coeff[j][k][(i+index)%m_filter_order] += place_holder*m_error[j]*m_x[k][(i+index)%m_filter_order]*(!m_triggered[k]);
					 if (i == m_step_size - 1 and j != k)
						 m_coeff[j][k][m_filter_order] += place_holder*m_error[j]*data[k]*(!m_triggered[k]);
				}
		}
	// Check for shot and increment variance
		for (int i = 0; i < m_buffer_number; i++)
		{
			m_triggered[i] = (m_error[i]*m_error[i] > m_threshold*m_var[i]);

			m_var[i] = m_alpha*m_var[i] + (1-m_alpha)*m_error[i]*m_error[i];
		// Add data to m_x for next round; If triggered, does not add that data point
			m_x[i][index%m_filter_order] = data[i]*m_triggered[i];
		}
		return true;
	}
*/
}

bool RingBuffer::IsTriggered(int &which_channel)
{
	bool val = false;
	for (int i = 0; i< m_buffer_number; i++)
	{
		if (m_triggered[i] == true)
		{
			which_channel = i;
			val = true;
			break;
		}
	}
	return val;
}

void RingBuffer::ResetTriggers()
{
	for(int i = 0; i < m_buffer_number; i++)
		m_triggered[i] = false;
}

bool RingBuffer::IsRecording()
{
	bool val = true;
	for (int i = 0; i< m_buffer_number; i++)
	{
		if (!m_record_history[i])
		{
			val = false;
			break;
		}
	}
	return val;
}
