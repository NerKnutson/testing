#ifndef WELCH_H
#define WELCH_H
#include <cmath>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>

class Welch
{
	public:
		Welch(size_t total_size, size_t window_size, size_t sampling_rate, size_t N_chan, double alpha) :
			m_total_size(total_size),
			m_alpha(alpha),
			m_N_chan(N_chan),
			m_window_size(window_size),
			m_window_num(2*total_size/window_size - 1),
			m_bin_num(window_size/2 - 1)
			{
				matrix_holder = new gsl_matrix_complex*[m_bin_num];
				for (int b = 0; b < m_bin_num; b++)
					matrix_holder[b] = gsl_matrix_complex_calloc(N_chan, N_chan);
				vector_holder = new gsl_vector_complex*[m_bin_num];
				for (int b = 0; b < m_bin_num; b++)
					vector_holder[b] = gsl_vector_complex_calloc(N_chan);
			}
		~Welch()
			{
				for(int b = 0; b < m_bin_num; b++)
				{
					gsl_matrix_complex_free(matrix_holder[b]);
					gsl_vector_complex_free(vector_holder[b]);
				}
				delete matrix_holder;
				delete vector_holder;
			}
		void tukey(double* input);
		void welch_csd(double* input);
		void fft(double* input);
	// Array of N_chan by N_chan sized matrices holding CSD for each frequency
		gsl_matrix_complex** matrix_holder;
	// Array of N_chan sized vectors holding FFT coefficients for each frequency
		gsl_vector_complex** vector_holder;
	private:
	// Total size of data set
		const size_t m_total_size;
	// Alpha for the windowing function
		const double m_alpha;
	// Number of channels
		const size_t m_N_chan;
	// Size of the window for the Welch method and size of the shot data set
		const size_t m_window_size;
	// Number of m_window_size windows contained in the m_total size history data set
		const size_t m_window_num;
	// Number of bins (number of FFT coefficients)
		const size_t m_bin_num;
};
#endif
