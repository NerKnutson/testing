#include "welch.h"

int main()
{
	size_t total_size = 1024;
	size_t window_size = 512;
	size_t N_chan = 1;
	double alpha = 0.5;
	Welch example(total_size,window_size,N_chan,alpha);

	double array[total_size] = {1};
	for (int i = 0; i < total_size; ++i)
		array[i] = 1.0;
	example.tukey((double*)array);
	printf("\n\n");
	for (int i = 0; i < window_size; ++i)
		printf("%10.5f\n",array[i]);
}
