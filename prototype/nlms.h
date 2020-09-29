#ifndef NLMS_H
#define NLMS_H

class NLMS
{
  private:
    double error;
    double* coeff;
    int m_filter_order;
    double m_step_size;
    double m_N_chan;

  public:

    NLMS(int filter_order, double step_size);
    double Process(const RingBuffer rbuffer);

};

#endif
