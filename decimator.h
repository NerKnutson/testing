#ifndef DECIMATOR_H
#define DECIMATOR_H
# include "Iir.h"
# include <cmath>
class Decimator
{
	private:
// Passband and stopband attenuations
		double m_gpass;
		double m_gstop;
// Passband and stopband edge frequencies
		double m_w_p;
		double m_w_s;
// Decimation Factor
		int m_df;
// Sampling Frequency
		double m_fs;
// Decimation state
// if (m_d_state == m_df) <DECIMATE DATA POINT>
		int m_d_state;
// Number of Channels
		int m_N_chan;
// Filter for each channel
		Iir::ChebyshevII::LowPass<10> *m_filt;
	public:
	// Public data
		double *m_data_out;
	// Parameterized constructor can create N channels
		Decimator(double gpass, double gstop, double w_p, double w_s, int d_factor, int N_chan, double fs_in):
			m_gpass(gpass),
			m_gstop(gstop),
			m_w_p(w_p),
			m_w_s(w_s),
			m_df(d_factor),
			m_fs(fs_in),
			m_d_state(0),
			m_N_chan(N_chan)
			{
				m_data_out = new double[m_N_chan];
				m_filt = new Iir::ChebyshevII::LowPass<10>[m_N_chan];
				for (int i = 0; i < N_chan; ++i) m_filt[i].setup(m_fs, m_w_p, m_gstop);
			}
	// Destructor
		~Decimator()
			{
				delete m_data_out;
				delete m_filt;
			}
	// Public methods
		bool ProcessData(double data[]);
		int GetOrder();
};
#endif
