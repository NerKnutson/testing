#ifndef DECIMATOR_H
#define DECIMATOR_H
# include "Iir.h"
# include <math.h>
class Decimator
{
	private:
		double m_gpass;
		double m_gstop;
		double m_w_p;
		double m_w_s;
		int m_df;
		double m_fs;
		int m_d_state;
		int m_N_chan;
		Iir::ChebyshevII::LowPass<10> *m_filt;
	public:
	// Public data
		double *m_data_out;
	// Parameterized constructor can create N channels
		Decimator(double gpass, double gstop, double w_p, double w_s, int d_factor, int N_chan, double fs_in);
	// Destructor
		~Decimator();
	// Public methods
		bool ProcessData(double data[]);
		int GetOrder();
};
#endif
