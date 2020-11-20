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

#include <chrono>


const size_t subarray_size = 4;

const size_t window_size = 256;
const size_t sampling_rate = 48000;
const double tukey_window_factor = 0.5;
const size_t N_chan = 19;

const double gpass = 0.1;
const double gstop = 40.0;
const double w_p = 1000.0;
const double w_s = 2000.0;
const size_t d_factor = 10;

const size_t history_size = 512;

const double threshold = 500;

void fourier_transformer(int threaded_channel, double* threaded_history, double* threaded_shot)
{
auto start = std::chrono::steady_clock::now();
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
//gsl_vector_fprintf(stdout,find_w_opt_params->slowness, "%10.5f");
//printf("\nSpeed:%10.5f\n\n", 1.0/gsl_blas_dnrm2(find_w_opt_params->slowness));

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> thread_time = end - start;
	printf("Time from start of thread: %10.5f\n",
			thread_time.count()
		);
}

int main()
{
//	Decimator things
	Decimator deci(gpass,gstop,w_p,w_s,d_factor,N_chan,sampling_rate);
	//printf("%d\n",deci.GetOrder());

	int which_channel = 0;
	//NOTE: we subtract the length of the part that would be distorted by the Tukey windowing function.
	size_t waiting_samples = window_size - floor(0.5*tukey_window_factor*(window_size+1));
	//size_t waiting_samples = window_size;
	int waiting_index = 0;
	bool start_waiting = false;

//	Ring Buffer
	RingBuffer rbuff(N_chan,history_size+window_size);

	std::ifstream infile;
	//infile.open("/home/ner/programming/inputs/signal_face_centered_cubic_x_axis.dat");
	//infile.open("/home/ner/programming/inputs/lower_noise_z_axis.dat");
	infile.open("/home/ner/programming/inputs/recording_1.dat");
	std::string line;
	std::string::size_type sz;
	double array[N_chan] = {0};

//	NLMS Variables
	const int filter_order = 3;
	const double step_size = 0.2;

	double alpha = 0.9999;
	const double tiny = 0.00000000001;

	double* array_setup[N_chan];

	FILE* a_g = fopen("array_groups","r");
	int temp = 0;
	for(int n = 0; n < N_chan; n++)
	{
		fscanf(a_g,"%d",&temp);
		array_setup[n] = &(deci.m_data_out[temp-1]);
	}

	NLMS filter_1(5, filter_order, step_size, alpha, threshold, tiny);
	NLMS filter_2(5, filter_order, step_size, alpha, threshold, tiny);
	NLMS filter_3(4, filter_order, step_size, alpha, threshold, tiny);
	NLMS filter_4(4, filter_order, step_size, alpha, threshold, tiny);
	NLMS filter_raw(1, filter_order, step_size, alpha, threshold, tiny);
	double* pointer;


int counter = 0;
int shot_num = 0;
if (infile.is_open())
{
	auto start = std::chrono::steady_clock::now();
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


			pointer = array_setup[0];
			filter_1.Process(pointer);

			pointer = array_setup[5];
			filter_2.Process(pointer);

			pointer = array_setup[10];
			filter_3.Process(pointer);

			pointer = array_setup[14];
			filter_4.Process(pointer);

			pointer = array_setup[18];
			filter_raw.Process(pointer);

			counter++;
		}

		//if (rbuff.Size() > history_size + floor(0.5*tukey_window_factor*(window_size+1)) and !start_waiting and n_filter.IsTriggered(which_channel))
		if (rbuff.Size() > history_size + floor(0.5*tukey_window_factor*(window_size+1))
				and !start_waiting
				and (filter_1.IsTriggered(which_channel)
					or filter_2.IsTriggered(which_channel)
					or filter_3.IsTriggered(which_channel)
					or filter_4.IsTriggered(which_channel)
					or filter_raw.IsTriggered(which_channel)
					)
				)
		//if (rbuff.Full() and !start_waiting and n_filter.IsTriggered(which_channel))
		{
			//printf("Data Point Triggered: %d\n",counter);
			start_waiting = true;

			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> trigger_time = end - start;
			printf("Expected time from start to trigger: %10.5f\n",
					1.0*counter/(sampling_rate/d_factor)
				);
			printf("Actual time from start to trigger: %10.5f\n\n",
					trigger_time.count()
				);

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
			filter_1.ResetTriggers();
			filter_2.ResetTriggers();
			filter_3.ResetTriggers();
			filter_4.ResetTriggers();
			filter_raw.ResetTriggers();

			rbuff.Dump(history_size,history,shot);
			rbuff.UnPut(window_size);
			std::thread (fourier_transformer,which_channel,history,shot).detach();
		}
	}
}
	infile.close();
	return 0;
}
