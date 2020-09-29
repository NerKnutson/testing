#include "rbuffer.h"
#include "nlms.h"
#include <iostream>

NLMS::NLMS(int filter_order, double step_size)
{
  m_filter_order = filter_order;
  m_step_size = step_size;

  error = 0;
  coeff = new double[m_filter_order];
}

double NLMS::Process(const RingBuffer rbuffer)
{
  error = rbuffer.Look(m_filter_order)[1];

  double place_holder;
  double tiny = 0.0000000001;


  for (int i = 0; i < m_filter_order; i++)
  {
    error -= coeff[i]*rbuffer.Look(i)[1];
    place_holder += rbuffer.Look(i)[1]*rbuffer.Look(i)[1];
  }

  place_holder += tiny;
  place_holder = m_step_size*error/place_holder;

  for (int i = 0; i < m_filter_order; i++)
  {
    std::cout << coeff[i] << "\n";
    coeff[i] = coeff[i] + place_holder;
    std::cout << coeff[i] << "\n";
  }

  return 0.0;
}
