#include <iostream>
#include <fstream>
#include "rbuffer.h"
#include <iomanip>
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


//  NLMS Variables
  const int m_filter_order = 7;
  const double m_step_size = 0.0001;

  double place_holder;
  const double tiny = 0.00000000001;

  double m_coeff[N_chan][N_chan][m_filter_order + 1] = {0};
  double m_error[N_chan] = {0};

  RingBuffer m_x(N_chan,m_filter_order);

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
      m_x.Put(array);

      if(m_x.Full())
      {
        place_holder = tiny;
        m_x.Look(m_filter_order,m_error);

        for (int i = 0; i <= m_filter_order; i++)
        {
          m_x.Look(i,array);
          for (int j = 0; j < N_chan; j++)
          {
              place_holder += array[j]*array[j];
              for (int k = 0; k < N_chan; k++)
              {
                if(i < m_filter_order)
                  m_error[j] -= array[k]*m_coeff[j][k][i];
                else if(j != k)
                  m_error[j] -= array[k]*m_coeff[j][k][i];
              }
          }
        }

        place_holder = m_step_size/place_holder;

        for (int i = 0; i <= m_filter_order; i++)
        {
          m_x.Look(i,array);
          for (int j = 0; j < N_chan; j++)
            for (int k = 0; k < N_chan; k++)
              m_coeff[j][k][i] += place_holder*m_error[j]*array[k];
        }

        for (int j = 0; j < N_chan; j++)
          std::cout << std::setw(10) << m_error[j];

        std::cout << std::endl;
      }
    }
  }



/*
for (int OVERALL = 0; OVERALL < rbuff.Size() - m_filter_order; OVERALL++)
{
  rbuff.Look(m_filter_order+OVERALL,error);

  for (int j = 0; j < N_chan; j++)
  {
    x[m_filter_order][j] = error[j];

    for (int k = 0; k < N_chan; k++)
    {
        if (k != j)
          error[k] -= error[j]*coeff[k][j][m_filter_order];
    }
  }

  place_holder = tiny;

  for (int i = 0; i < m_filter_order; i++)
  {
    rbuff.Look(i+OVERALL,x[i]);

    for (int j = 0; j < N_chan; j++)
    {
      for (int k = 0; k < N_chan; k++)
      {
        error[k] -= x[i][j]*coeff[k][j][i];
      }
      place_holder += x[i][j]*x[i][j];
    }
  }

  place_holder = m_step_size/place_holder;

  for(int j = 0; j < N_chan; j++)
    for(int k = 0; k < N_chan; k++)
    {
      for(int i = 0; i <= m_filter_order; i++)
      {
        coeff[k][j][i] += place_holder*error[k]*x[i][j];
      }
    }

  for (int k = 0; k < N_chan; k++)
    std::cout << std::setw(10) << error[k];

  std::cout << std::endl;
  break;
}
*/




/*
for (int i = 0; i < m_filter_order; i++)
{
  rbuff.Look(i,x[i]);
  if (i < N_chan)
    error[i] = x[0][i];

  for (int j = 0; j < N_chan && i != 0; j++)
  {
    for (int k = 0; k < N_chan; k++)
      error[k] -= x[i][j]*coeff[j][k*m_filter_order + i];
  }
}

for (int j = 0; j < N_chan; j++)
{
  for (int i = 0; i < m_filter_order; i++)
    std::cout << std::setw(10) << x[i][j]<< "\t";

  std::cout << std::endl;
}
*/


}
