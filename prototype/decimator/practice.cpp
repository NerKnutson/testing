#include <iostream>
#include <fstream>
#include <string>
#include "decimator.h"
#include "rbuffer.h"

int main()
{

  // decimator parameters
  double w_p = 1200.0;
  double w_s = 1600.0;
  double gpass = 0.1;
  double gstop = 40.0;
  int d_factor = 10;
  int N_chan = 5;
  double fs_in = 48000.0;

  // create decimator
  Decimator deci(gpass, gstop, w_p, w_s, d_factor, N_chan, fs_in);

  // Create Ring Buffer
  RingBuffer rbuff(N_chan,10000);

  std::ifstream infile;
  infile.open("input.dat");
  std::string line;
  std::string::size_type sz;
  double array[N_chan] = {0};

  // loop over data file
  if (infile.is_open())
    {
      while (getline(infile, line))
	{
	  // delete first and last characters
	  //line.erase(line.begin());
	  //line.erase(line.end() - 1);

	  // load an array to feed the decimator
	  for (int i = 0; i < N_chan; ++i)
	    {
	      array[i] = std::stod(line, &sz);
	      line = line.substr(sz);
	    }

	  // call the decimator, and if it is ready, print to stdout
	  if (deci.ProcessData(array))
	    {
         // insert decimator output into ring buffer
         rbuff.Put(deci.m_data_out);

/*
         // retrieve contents of ring buffer and print
         double* output;
         output=rbuff.Get();
         for(int i=0;i<N_chan;i++)
         {
            std::cout << output[i] << " ";
         }
*/
	}
    }


  double*** output;
  output = rbuff.Dump(100);

  for(int i=0;i<100;i++)
  {
    for(int j=0;j<N_chan;j++)
    {
      std::cout << output[0][i][j] << "\t";
    }
    std::cout << "\n";
  }

  for(int i=100;i<10000;i++)
  {
    for(int j=0;j<N_chan;j++)
    {
      std::cout << output[1][i][j] << "\t";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  }

  // Print order of decimator
  std::cout << "\n" << deci.GetOrder() << std::endl;
}
