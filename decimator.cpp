#include "decimator.h"

// Parameterized constructor can create N channels
Decimator::Decimator(double gpass, double gstop, double w_p, double w_s, int d_factor, int N_chan, double fs_in)
{
	m_gpass = gpass;
	m_gstop = gstop;
	m_w_p = w_p;
	m_w_s = w_s;
	m_df = d_factor;
	m_fs =  fs_in;
	m_d_state = 0;
	m_N_chan = N_chan;
	m_data_out = new double[m_N_chan];
	m_filt = new Iir::ChebyshevII::LowPass<10>[m_N_chan];
	for (int i = 0; i < N_chan; ++i) m_filt[i].setup(m_fs, m_w_p, m_gstop);
}

// Destructor
Decimator::~Decimator()
{
	delete m_filt;
	delete m_data_out;
}

// Decimator
bool Decimator::ProcessData(double data[])
{
	for(int i = 0; i < m_N_chan; ++i)
	 {
		m_data_out[i] = m_filt[i].filter(data[i]);
	 }

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
