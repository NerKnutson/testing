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

const size_t window_size = 256;
const size_t sampling_rate = 48000;
const double tukey_window_factor = 0.25;
const size_t N_chan = 19;

const double gpass = 0.1;
const double gstop = 40.0;
const double w_p = 1000.0;
const double w_s = 2000.0;
const size_t d_factor = 10;

const size_t history_size = 1024;

const double threshold = 500;

void fourier_transformer(int threaded_channel, double* threaded_history, double* threaded_shot)
{
// Creates Fourier object and processes history and shot data
	Welch fourier_stuff(history_size,window_size,N_chan,tukey_window_factor);
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
//ntk_test_routine(find_w_opt_params);
gsl_vector_fprintf(stdout,find_w_opt_params->slowness, "%10.5f");
fflush(stdout);
printf("\nSpeed:%10.5f\n\n", 1.0/gsl_blas_dnrm2(find_w_opt_params->slowness));
fflush(stdout);

}

int main()
{
//	Decimator things
	Decimator deci(gpass,gstop,w_p,w_s,d_factor,N_chan,sampling_rate);
	//printf("%d\n",deci.GetOrder());
//	Waiting for after trigger before dumping
	double t_c = t_d + t_t;	//	(ms)

	int which_channel = 0;
	//size_t waiting_samples = t_c*sampling_rate/1000;
	//NOTE: we subtract the length of the part that would be distorted by the Tukey windowing function.
	int tukey_length = floor(0.5*tukey_window_factor*(window_size+1));
	size_t waiting_samples = window_size - tukey_length;
	//size_t waiting_samples = window_size;
	int waiting_index = 0;
	bool start_waiting = false;

//	Ring Buffer
	RingBuffer buff(N_chan,history_size + window_size);

	std::ifstream infile;
	//infile.open("/home/ner/programming/inputs/signal_face_centered_cubic_x_axis.dat");
	//infile.open("/home/ner/programming/inputs/lower_noise_z_axis.dat");
	//infile.open("/home/ner/programming/inputs/recording_1.dat");
	infile.open("/home/ner/programming/inputs/synthetic_zylia.dat");
	std::string line;
	std::string::size_type sz;
	double array[N_chan] = {0};

//	NLMS Variables
	const int filter_order = 10;
	const double step_size = 0.2;

	double alpha = 0.9999;
	const double tiny = 0.00000000001;

	NLMS n_filter(N_chan, filter_order, step_size, alpha, threshold, tiny);


int counter = 0;
int shot_num = 0;
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
			n_filter.Process((double*)deci.m_data_out);
			buff.Put((double*)deci.m_data_out);
			counter++;
		}

		if (buff.Size() > tukey_length and !start_waiting and n_filter.IsTriggered(which_channel))
		{
			printf("Data Point Triggered: %d\n",counter);
			fflush(stdout);
			start_waiting = true;
		}

		if (start_waiting == true)
			waiting_index++;

		if (waiting_index/d_factor == waiting_samples)
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

			buff.Dump(history_size,history,shot);
			buff.UnPut(window_size);
			std::thread (fourier_transformer,which_channel,history,shot).detach();
/*
			for(int w = 0; w < window_size; w++)
			{
				for(int n = 0; n < N_chan; n++)
					printf("%20.15f", shot[n*window_size + w]);
				printf("\n");
			}
			printf("\n\n");
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
	return 0;
}
