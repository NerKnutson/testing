#include <iostream>
#include <fstream>
#include "rbuffer.h"
//#include "nlms.h"

int main()
{
  int N_chan = 5;
  RingBuffer rbuff(N_chan,1000000);

  std::ifstream infile;
  infile.open("only_sine.dat");
  std::string line;
  std::string::size_type sz;
  double array[N_chan] = {0};

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
    }
  }

  int m_filter_order = 10;
  double m_step_size = 0.0001;

  double place_holder;
  double tiny = 0.0000000001;

  double coeff[m_filter_order] = {0};
  double error = 0;
  double x[N_chan];

  int channel = 1;

for (int j = 0; j < rbuff.Size()-m_filter_order; j++)
{
  rbuff.Look(m_filter_order+j,x);
  error = x[channel];

  place_holder = tiny;

  for (int i = 0; i < m_filter_order; i++)
  {
    rbuff.Look(i+j,x);
    error -= coeff[i]*x[channel];
    place_holder += x[channel]*x[channel];
  }

  place_holder = m_step_size*error/place_holder;

  for (int i = 0; i < m_filter_order; i++)
  {
    rbuff.Look(i+j,x);
    coeff[i] = coeff[i] + place_holder*x[channel];
    //std::cout << coeff[i] << "\t";
  }

//  cout << "\n\nActual:\t" << x[channel]
//       << "\tError:\t" << error << "\n\n";
  std::cout << error << std::endl;
}



  return 0;
}
