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
		//  ERROR = PRESENT_DATA
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
		place_holder = m_step_size/place_holder;
		for (int i = 0; i < m_filter_order; i++)
		{
			for (int j = 0; j < m_N_chan; j++)
				for (int k = 0; k < m_N_chan; k++)
				{
					 m_coeff[j][k][(i+index)%m_filter_order] += place_holder*m_error[j]*m_x[k][(i+index)%m_filter_order];
					 if (i == m_step_size - 1)
						 m_coeff[j][k][m_filter_order] += place_holder*m_error[j]*data[k];
				}
		}
	// Check for shot and increment variance
		for (int i = 0; i < m_N_chan; i++)
		{
			if ( m_error[i]*m_error[i] > 25*m_var[i])
				m_triggered[i] = true;
			else
				m_triggered[i] = false;
			m_var[i] = m_alpha*m_var[i] + (1-m_alpha)*m_error[i]*m_error[i];
		// Add data to m_x for next round
			m_x[i][index%m_filter_order] = data[i];
		}
		return true;
	}
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
