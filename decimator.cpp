#include "decimator.h"

// Decimator
bool Decimator::ProcessData(double data[])
{
// Low Pass Filters the datapoint in
	for(int i = 0; i < m_N_chan; ++i)
	 {
		m_data_out[i] = m_filt[i].filter(data[i]);
	 }

// Tests if that point is decimated or not
	m_d_state += 1;

	if (m_d_state == m_df)
	 {
		m_d_state = 0;
		return true;
	 }
	else
		return false;
	}

// Calculate and Return Order
int Decimator::GetOrder()
{
	return ceil(acosh(sqrt( (pow(10,0.1*m_gstop)-1)/(pow(10,0.1*m_gpass)-1) ))/acosh(m_w_s/m_w_p));
}
