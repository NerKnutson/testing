#include "nlms.h"

bool NLMS::Process(double data[])
{
	if (index < m_filter_order)
	{
	// Fill array with data until we have enough to begin NLMS
		for (int i = 0; i < m_N_chan; i++)
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
		for (int i = 0; i < m_N_chan; i++)
		{
			m_error[i] = data[i];
			place_holder += data[i]*data[i];
		}
	// ERROR -= OLD_DATA
		for (int i = 0; i < m_filter_order; i++)
		{
			for (int j = 0; j < m_N_chan; j++)
			{
				place_holder += m_x[j][i]*m_x[j][i];
				for (int k = 0; k < m_N_chan; k++)
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
			for (int j = 0; j < m_N_chan; j++)
				for (int k = 0; k < m_N_chan; k++)
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
		for (int i = 0; i < m_N_chan; i++)
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
		for (int i = 0; i < m_N_chan; i++)
		{
			m_error[i] = data[i]*(!m_triggered[i]);
			place_holder += data[i]*data[i]*(!m_triggered[i]);
		}
	// ERROR -= OLD_DATA
		for (int i = 0; i < m_filter_order; i++)
		{
			for (int j = 0; j < m_N_chan; j++)
			{
				place_holder += m_x[j][i]*m_x[j][i]*(!m_triggered[j]);
				for (int k = 0; k < m_N_chan; k++)
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
			for (int j = 0; j < m_N_chan; j++)
				for (int k = 0; k < m_N_chan; k++)
				{
					 m_coeff[j][k][(i+index)%m_filter_order] += place_holder*m_error[j]*m_x[k][(i+index)%m_filter_order]*(!m_triggered[k]);
					 if (i == m_step_size - 1 and j != k)
						 m_coeff[j][k][m_filter_order] += place_holder*m_error[j]*data[k]*(!m_triggered[k]);
				}
		}
	// Check for shot and increment variance
		for (int i = 0; i < m_N_chan; i++)
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

bool NLMS::IsTriggered(int &which_channel)
{
	bool val = false;
	for (int i = 0; i< m_N_chan; i++)
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

void NLMS::ResetTriggers()
{
	for(int i = 0; i < m_N_chan; i++)
		m_triggered[i] = false;
}

bool NLMS::IsRecording()
{
	bool val = true;
	for (int i = 0; i< m_N_chan; i++)
	{
		if (!m_record_history[i])
		{
			val = false;
			break;
		}
	}
	return val;
}
