#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>

const double alpha = 0.75;

const size_t total_size = 4096;
const size_t N_chan = 5;
const size_t window_size = 256;
const size_t window_num = 2*total_size/window_size - 1;

const size_t bin_num = 4096*1000/total_size;
gsl_matrix_complex* matrix_holder[bin_num];


void tukey(double* input)
{
	double filter = 0.0;
	for (int i = 0; i < 0.5*alpha*window_size; i++)
	{
		filter = 0.5*(1-cos(2*M_PI*i/(alpha*window_size)));
		input[i] = input[i]*filter;
		input[window_size - i] = input[window_size - i]*filter;
	}
}

void fft(double* input)
{
	double psd = 0.0;
	double data[N_chan][window_size];
	size_t window_num = 2*total_size/window_size - 1;
	for (int b = 0; b < window_num; b++)
	{
		for (int n = 0; n < N_chan; n++)
		{
			for (int i = 0; i < window_size; i++)
			{
				data[n][i] = input[2*b*window_size + i*N_chan + n];
			}
			tukey(data[n]);
			gsl_fft_real_radix2_transform(data[n],1,window_size);
		}
		for (int w = 0; w < bin_num; w++)
		for (int n = 0; n < N_chan; n++)
			for (int m = n; m < N_chan; m++)
			{
				gsl_matrix_complex_set(
					matrix_holder[w], n, m,
					gsl_complex_add(
						gsl_matrix_complex_get(matrix_holder[0],n,m),
						gsl_complex_rect(
							data[n][w]*data[m][w] + data[n][window_size-1-w]*data[m][window_size-1-w],
							-data[n][w]*data[m][window_size-1-w] + data[n][window_size-1-w]*data[m][w]
						)
					)
				);
				if (b == window_num-1)
					gsl_matrix_complex_set(
						matrix_holder[w],n,m,
						gsl_complex_div_real(gsl_matrix_complex_get(matrix_holder[w],n,m),window_num)
					);
			}
	}
}


int main ()
{
	std::ifstream infile;
	infile.open("input.dat");
	std::string line;
	std::string::size_type sz;

	double input[total_size][N_chan];
	for (int w = 0; w < bin_num; w++)
		matrix_holder[w] = gsl_matrix_complex_calloc(N_chan, N_chan);

	size_t i = 0;
	if (infile.is_open())
		while (getline(infile,line))
		{
			for (int n = 0; n < N_chan; n++)
			{
				input[i][n] = std::stod(line,&sz);
				line = line.substr(sz);
			}
			i++;
		}
	fft((double*)input);
	for (int w = 0; w < bin_num; w++)
	{
		for (int n = 0; n < N_chan; n++)
		{
			std::cout << std::setw(16) << GSL_REAL(gsl_matrix_complex_get(matrix_holder[w],n,n));
		}
		std::cout << std::endl;
	}

/*
	for (int w = 0; w < bin_num; w++)
	{
		for (int n = 0; n < N_chan; n++)
		{
			for (int m = 0; m < N_chan; m++)
			{
				std::cout << std::setw(16) << GSL_IMAG(gsl_matrix_complex_get(matrix_holder[w],n,m));
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
*/
	infile.close();
	for (int w = 0; w < bin_num; w++)
		gsl_matrix_complex_free(matrix_holder[w]);
}
