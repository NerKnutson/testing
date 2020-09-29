#include <iostream>     //  cout, endl
#include <fstream>      //  ifstream
#include <iomanip>      //  setw
#include "rbuffer.h"    // RingBuffer

int main()
{
  int N_chan = 5;
  size_t buffer_size = 100000;

  RingBuffer rbuffer(N_chan, buffer_size);

  std::ifstream infile;
  infile.open("fuzzy_sine.dat");
  std::string line;
  std::string::size_type sz;
  double array[N_chan] = {0};

//  Filling RingBuffer
  if (infile.is_open())
  {
    while (getline(infile, line))
    {
      for (int i = 0; i < N_chan; i++)
      {
        array[i] = std::stod(line, &sz);
        line = line.substr(sz);
      }
      rbuffer.Put(array);

      if(rbuffer.Full())
        break;
    }
  }
  std::cout << "Putting is done" << std::endl;

//  Looking at RingBuffer
  double looking[N_chan];

  for (int i = 0; i < rbuffer.Size(); i++)
  {
    rbuffer.Look(i,looking);

    for (int j = 0; j < N_chan; j++)
      std::cout << std::setw(15) << looking[j];

    std::cout << std::endl;
  }

  std::cout << "Looking is done" << std::endl;

//  Dumping RingBuffer
  int marker = 1000;
  double* history[marker];
  double* shot[buffer_size-marker];

  //  Constructing arrays to Dump into
  for (int i = 0; i < marker; i++)
    history[i] = new double[N_chan];

  for (int i = marker; i < buffer_size; i++)
    shot[i-marker] = new double[N_chan];

  rbuffer.Dump(marker, history, shot);

  //  Showing what was dumped
  for (int i = 0; i < marker; i++)
  {
    for (int j = 0; j < N_chan; j++)
      std:: cout << std::setw(15) << history[i][j];

    std::cout << std::endl;
  }

  for (int i = marker; i < buffer_size; i++)
  {
    for (int j = 0; j < N_chan; j++)
      std:: cout << std::setw(15) << shot[i-marker][j];

    std::cout << std::endl;
  }

  std::cout << "Dumping is done" << std::endl;

//  Getting data from buffer

  double getting[N_chan];

  while (!rbuffer.Empty())
  {
    rbuffer.Get(getting);
    for (int i = 0; i < N_chan; i++)
      std::cout << std::setw(15) <<  getting[i];

    std::cout << std::endl;
  }

  std::cout << "Getting is done" << std::endl;

}
