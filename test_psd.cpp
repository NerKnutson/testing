#include <iostream>
#include <fstream>
#include <iomanip>
#include "rbuffer.h"
#include "nlms.h"
#include "welch.h"
#include "decimator.h"

#include <gsl/gsl_fft_complex.h>

#define t_d 75	//	length of signal in (ms)
#define t_t 25	//	transient time to cross array in (ms)

const size_t total_size = 8192;
const size_t window_size = 256;
const size_t sampling_rate = 10000;
const double alpha_2 = 0.50;
const size_t N_chan = 5;

const double gpass = 1.0;
const double gstop = 1.0;
const double w_p = 80;
const double w_s =100;
const size_t d_factor = 10;

int main()
{
//	Ring Buffer
	RingBuffer rbuff(N_chan,total_size);

	std::ifstream infile;
	infile.open("data/noise.dat");
	std::string line;
	std::string::size_type sz;
	double array[N_chan] = {0};

//	Welch Things
	double dumped[N_chan][total_size];
	Welch foo(total_size,window_size,sampling_rate,N_chan,alpha_2);

int linenumber = 0;
	if (infile.is_open())
	{
		while (getline(infile, line))
		{
			for (int i = 0; i < N_chan; i++)
			{
			array[i] = std::stod(line, &sz);
			line = line.substr(sz);
			}

				linenumber++;
				rbuff.Put(array);
			if (rbuff.Full())
				break;
		}
	}
	rbuff.Dump(total_size,(double*)dumped,(double*)array);
	foo.fft((double*)dumped);
	for (int b = 0; b < window_size/2 - 1; b++)
	{
		for (int n = 0; n < N_chan; n++)
		{
			std::cout << std::setw(16) << GSL_REAL(gsl_matrix_complex_get(foo.matrix_holder[b],n,n));
/*
			for (int m = 0; m < N_chan; m++)
			{
				std::cout << std::setw(16) << GSL_REAL(gsl_matrix_complex_get(foo.matrix_holder[b],n,m));
			}
			std::cout << std::endl;
*/
		}
		std::cout << std::endl;
	}
}
