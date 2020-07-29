#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <gsl/gsl_fft_real.h>

const double alpha = 0.75;

const size_t total_size = 4096;
const size_t N_chan = 5;
const size_t window_size = 256;


void tukey(size_t N_chan, size_t window_size, double* data)
{
	double filter = 0.0;
	for (int i = 0; i < 0.5*alpha*window_size; i++)
	{
		for (int n = 0; n < N_chan; n++)
		{
			filter = 0.5*(1-cos(2*M_PI*i/(alpha*window_size)));
			data[n*window_size + i] = data[n*window_size + i]*filter;
			data[(n+1)*window_size - i] = data[(n+1)*window_size - i]*filter;
		}
	}
}
void split(size_t total_size, size_t window_size, size_t N_chan, const double* input, double* output)
{
	size_t bins = 2*total_size/window_size - 1;
	for (int b = 0; b < bins; b++)
	{
		for (int n = 0; n < N_chan; n++)
		{
			for (int i = 0; i < window_size; i++)
			{
				output[b*N_chan*window_size + n*window_size + i] = input[2*b*window_size + i*N_chan + n];
			}
		}
	}
}


int main ()
{
	std::ifstream infile;
	infile.open("random.dat");
	std::string line;
	std::string::size_type sz;

	size_t bins = 2*total_size/window_size - 1;
	double data[total_size][N_chan];
	double split_data[bins][N_chan][window_size];

	size_t i = 0;
	if (infile.is_open())
		while (getline(infile,line))
		{
			for (int n = 0; n < N_chan; n++)
			{
				data[i][n] = std::stod(line,&sz);
				line = line.substr(sz);
			}
			i++;
		}
	split(total_size,window_size,N_chan,(double*)data,(double*)split_data);
	for (int b = 0; b < bins; b++)
	{
		tukey(N_chan,window_size,(double*)split_data[b]);
		for (int n = 0; n < N_chan; n++)
			gsl_fft_real_radix2_transform(split_data[b][n],1,window_size);
	}
	/*
		for (int i = 0; i < window_size; i++)
		{
			for (int n = 0; n < N_chan; n++)
				std::cout << std::setw(15) << data[n][i];
			std::cout << std::endl;
		}
	*/
	for (int b = 0; b < bins; b++)
	{
		for (int i = window_size; i > window_size/2; i--)
		{
			for (int n = 0; n < N_chan; n++)
				std::cout << std::setw(15) << split_data[b][n][i];
			std::cout << std::endl;
		}
	}
	infile.close();
}
