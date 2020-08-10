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
//	Decimator things
	Decimator deci(gpass,gstop,w_p,w_s,d_factor,N_chan,sampling_rate);
//	Waiting for after trigger before dumping
	double t_c = t_d + t_t;	//	(ms)

	int which_channel = 0;
	size_t waiting_samples = t_c*sampling_rate/1000;
	int waiting_index = 0;
	bool start_waiting = false;

//	Ring Buffer
	RingBuffer rbuff(N_chan,8192+waiting_samples);

	std::ifstream infile;
	infile.open("data/noise.dat");
	std::string line;
	std::string::size_type sz;
	double array[N_chan] = {0};

//	NLMS Variables
	const int filter_order = 7;
	const double step_size = 0.01;

	double alpha = 0.9999;
	const double tiny = 0.00000000001;

	NLMS n_filter(N_chan, filter_order, step_size, alpha, tiny);

//	Arrays made for dumping
	size_t history_size = rbuff.Capacity() - waiting_samples;
	double history[N_chan][history_size];
	double shot[N_chan][waiting_samples];

//	Welch Things
	Welch foo(history_size,window_size,sampling_rate,N_chan,alpha_2);

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
				rbuff.Put((double*)array);
				n_filter.Process((double*)array);
/*
			if (deci.ProcessData(array))
			{
				linenumber++;
				rbuff.Put((double*)deci.m_data_out);
				n_filter.Process((double*)deci.m_data_out);
			}
*/

			//if (rbuff.Full() and n_filter.IsTriggered(which_channel))
			if (rbuff.Full() and (linenumber==8193 or n_filter.IsTriggered(which_channel)))
			{
				start_waiting = true;
			}

			if (start_waiting == true)
				waiting_index++;

			if (waiting_index == waiting_samples)
			{
				rbuff.Dump(history_size,(double*)history,(double*)shot);
				foo.fft((double*)history);
/*
				for (int i = 0; i < history_size; i++)
				{
					for(int j = 0; j < N_chan; j++)
						std::cout << std::setw(15) << history[j][i];
					for(int j = 0; j < N_chan and i < waiting_samples; j++)
						std::cout << std::setw(15) << shot[j][i];
					std::cout << std::endl;
				}
*/
			}
		}
	}
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
