#include "welch.h"

void Welch::tukey(double* input)
{
	double filter = 0.0;
	for (int i = 0; i < 0.5*m_alpha*(m_window_size+1); i++)
	{
		filter = 0.5*(1-cos(2*M_PI*i/(m_alpha*(m_window_size+1))));
		input[i] = input[i]*filter;
		input[m_window_size-1 - i] = input[m_window_size-1 - i]*filter;
	}
}

void Welch::fft(double* input)
{
	double data[m_N_chan][m_window_size];
	for (int w = 0; w < m_window_num; w++)
	{
		for (int n = 0; n < m_N_chan; n++)
		{
			for (int i = 0; i < m_window_size; i++)
			{
				data[n][i] = input[w*m_window_size/2 + n*m_total_size + i];
			}
			tukey(data[n]);
			gsl_fft_real_radix2_transform(data[n],1,m_window_size);
		}
		for (int b = 0; b < m_bin_num; b++)
			for (int n = 0; n < m_N_chan; n++)
			for (int m = n; m < m_N_chan; m++)
			{
				gsl_matrix_complex_set(
					matrix_holder[b], n, m,
					gsl_complex_add(
						gsl_matrix_complex_get(matrix_holder[b],n,m),
						gsl_complex_rect(
							data[n][b+1]*data[m][b+1] + data[n][m_window_size-b-1]*data[m][m_window_size-b-1],
							-data[n][b+1]*data[m][m_window_size-b-1] + data[n][m_window_size-b-1]*data[m][b+1]
						)
					)
				);
				if (w == m_window_num-1)
				{
					gsl_matrix_complex_set(
						matrix_holder[b],n,m,
						gsl_complex_div_real(
							gsl_matrix_complex_get(
								matrix_holder[b],n,m
							),
							m_window_num
						)
					);
				}
		}
	}
}