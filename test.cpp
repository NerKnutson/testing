#include <iostream>
#include <fstream>
#include <iomanip>
#include <thread>
#include <mutex>

#include "rbuffer.h"
#include "nlms.h"
#include "welch.h"
#include "decimator.h"
#include "find_w.h"

#include <gsl/gsl_fft_complex.h>


#define t_d 75	//	length of signal in (ms)
#define t_t 25	//	transient time to cross array in (ms)

const size_t window_size = 1024;
const size_t sampling_rate = 48000;
const double tukey_window_factor = 0.25;
const size_t N_chan = 5;

const double gpass = 0.1;
const double gstop = 1;
const double w_p = 150.0;
const double w_s = 3000.0;
const size_t d_factor = 10;

const size_t history_size = 8192;

//void fourier_transformer(Welch threaded_fourier_stuff, double* threaded_history, double* threaded_shot)
void fourier_transformer(int threaded_channel, double* threaded_history, double* threaded_shot)
{
// Creates Fourier object and processes history and shot data
	Welch fourier_stuff(history_size,window_size,sampling_rate,N_chan,tukey_window_factor);
	fourier_stuff.welch_csd(threaded_history);
	fourier_stuff.fft(threaded_shot);

opt_params* find_w_opt_params;
int error = init_find_w("find_w.ini",&find_w_opt_params);
	if (error != 0)
	{
		printf("Initialization error: %d\n", error);
		exit(error);
	}

perform_opt(find_w_opt_params,fourier_stuff.matrix_holder,fourier_stuff.vector_holder,threaded_channel);

	printf("Finished optimization\n");
	fflush(stdout);
	for (int i = 0; i < 3; i++)
	{
			printf("%f \n", gsl_vector_get(find_w_opt_params->slowness,i) );
			fflush(stdout);
	}
}

int main()
{
//	Decimator things
	Decimator deci(gpass,gstop,w_p,w_s,d_factor,N_chan,sampling_rate);
//	Waiting for after trigger before dumping
	double t_c = t_d + t_t;	//	(ms)

	int which_channel = 0;
	//size_t waiting_samples = t_c*sampling_rate/1000;
	//NOTE: we subtract the length of the part that would be distorted by the Tukey windowing function.
	size_t waiting_samples = window_size - floor(0.5*tukey_window_factor*(window_size+1));
	int waiting_index = 0;
	bool start_waiting = false;

//	Ring Buffer
	RingBuffer rbuff(N_chan,history_size+window_size);

	std::ifstream infile;
	infile.open("input/recording_1.dat");
	std::string line;
	std::string::size_type sz;
	double array[N_chan] = {0};

//	NLMS Variables
	const int filter_order = 7;
	const double step_size = 0.2;

	double alpha = 0.9999;
	const double tiny = 0.00000000001;

	NLMS n_filter(N_chan, filter_order, step_size, alpha, tiny);

//	Arrays made for dumping
	//double history[N_chan][history_size];
	//double shot[N_chan][window_size];

//	Welch Things
	//Welch fourier_stuff(history_size,window_size,sampling_rate,N_chan,tukey_window_factor);

int counter = 0;

if (infile.is_open())
{
	while (getline(infile, line))
	{
		for (int i = 0; i < N_chan; i++)
		{
			array[i] = std::stod(line, &sz);
			line = line.substr(sz);
		}

		if (deci.ProcessData(array))
		{
			rbuff.Put((double*)deci.m_data_out);
			n_filter.Process((double*)deci.m_data_out);
			counter ++;
		}
		if (rbuff.Size() > history_size + floor(0.5*tukey_window_factor*(window_size+1)) and !start_waiting and n_filter.IsTriggered(which_channel))
		{
			start_waiting = true;
		}

		if (start_waiting == true)
			waiting_index++;

		if (waiting_index == waiting_samples)
		{
		// Arrays to dump in are made
			double* history;
				history = new double[N_chan*history_size];
			double* shot;
				shot = new double[N_chan*window_size];

			// Main thread dumps and resets triggers to begin collecting again
			waiting_index = 0;
			start_waiting = false;
			n_filter.ResetTriggers();
			rbuff.Dump(history_size,history,shot);

			printf("Starting new thread\n");
			fflush(stdout);
			std::thread (fourier_transformer,which_channel,history,shot).detach();
/*
			for(int w = 0; w < window_size; w++)
			{
				for(int n = 0; n < N_chan; n++)
					std::cout << std::setw(16) << shot[n][w];
				std::cout << std::endl;
			}
			for(int h = 0; h < history_size; h++)
			{
				for(int n = 0; n < N_chan; n++)
					std::cout << std::setw(16) << history[n][h];
				std::cout << std::endl;
			}
*/
		}
	}
}
	infile.close();
}
