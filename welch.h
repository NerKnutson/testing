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
			}
		~Welch()
			{
				for(int b = 0; b < m_bin_num; b++)
					gsl_matrix_complex_free(matrix_holder[b]);
				delete matrix_holder;
			}
		void fft(double* input);
		void tukey(double* input);
		gsl_matrix_complex** matrix_holder;
	private:
		const size_t m_total_size;
		const double m_alpha;
		const size_t m_N_chan;
		const size_t m_window_size;
		const size_t m_window_num;
		const size_t m_bin_num;
};
#endif
