#include <iostream>
#include <fstream>
#include <iomanip>
#include "rbuffer.h"
#include "nlms.h"

#include <gsl/gsl_fft_complex.h>

#define t_d 75	//	length of signal in (ms)
#define t_t 25	//	transient time to cross array in (ms)

int main()
{
	int N_chan = 5;
	RingBuffer rbuff(N_chan,4000);

	std::ifstream infile;
	infile.open("middle_spike.dat");
	std::string line;
	std::string::size_type sz;
	double array[N_chan] = {0};

//	NLMS Variables
	const int filter_order = 7;
	const double step_size = 0.01;

	double alpha = 0.9999;
	const double tiny = 0.00000000001;

	NLMS n_filter(N_chan, filter_order, step_size, alpha, tiny);

//	Waiting for after trigger before dumping
	double t_c = t_d + t_t;	//	(ms)
	double sampling_rate = 10000;	//	(Hz)

	int which_channel = 0;
	size_t waiting_samples = t_c*sampling_rate/1000;
	int waiting_index = 0;
	bool start_waiting = false;


//	Arrays made for dumping
	size_t history_size = rbuff.Capacity() - waiting_samples;
	double history[history_size][N_chan];
	double shot[waiting_samples][N_chan];


int linenumber = 0;
	if (infile.is_open())
	{
		while (getline(infile, line))
		{
			linenumber++;
			for (int i = 0; i < N_chan; i++)
			{
			array[i] = std::stod(line, &sz);
			line = line.substr(sz);
			}

			rbuff.Put(array);

			n_filter.Process(array);

			if (rbuff.Full() and n_filter.IsTriggered(which_channel))
				start_waiting = true;

			if (start_waiting == true)
				waiting_index++;

			if (waiting_index == waiting_samples)
			{
				rbuff.Dump(history_size,(double*)history,(double*)shot);
				for (int i = 0; i < history_size; i++)
				{
					for(int j = 0; j < N_chan; j++)
						std::cout << std::setw(15) << history[i][j];
					for(int j = 0; j < N_chan and i < waiting_samples; j++)
						std::cout << std::setw(15) << shot[i][j];
					std::cout << std::endl;
				}
			}
		}
	}
}
