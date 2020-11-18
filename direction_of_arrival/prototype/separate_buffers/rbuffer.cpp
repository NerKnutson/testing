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
