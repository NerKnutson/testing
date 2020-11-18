#include <iostream>
#include <fstream>
#include "rbuffer.h"
#include <iomanip>
#include "nlms.h"

int main()
{
  int N_chan = 5;
  RingBuffer rbuff(N_chan,10000);

  std::ifstream infile;
  infile.open("middle_spike.dat");
  std::string line;
  std::string::size_type sz;
  double array[N_chan] = {0};


//  NLMS Variables
  const int filter_order = 7;
  const double step_size = 0.01;

  double alpha = 0.9999;
  const double tiny = 0.00000000001;

  NLMS n_filter(N_chan, filter_order, step_size, alpha, tiny);


int entry_number = 0;

  if (infile.is_open())
  {
    while (getline(infile, line))
    {
      for (int i = 0; i < N_chan; i++)
      {
        array[i] = std::stod(line, &sz);
        line = line.substr(sz);
      }
      rbuff.Put(array);


      if(n_filter.Process(array))
      {
        for (int j = 0; j < N_chan; j++)
        {
          if (n_filter.m_triggered[j])
          {
            if(j==0)
              std::cout << std::setw(15) << entry_number << std::setw(15) << n_filter.m_var[j];
            else if(j!= N_chan-1)
              std::cout << std::setw(15) << n_filter.m_var[j];
            else
              std::cout << std::setw(15) << n_filter.m_var[j] << std::endl;
          }
        }
      }
      else
        std::cout << "Not enough data" << std::endl;
    }
  }
}
