#ifndef RBUFFER_H
#define RBUFFER_H
#include <mutex>
#include <memory>
#include <cstddef>

class RingBuffer
{
	public:
		explicit RingBuffer(size_t buffer_number, size_t size):
			m_buffer_number(buffer_number),
			max_size_(size)
			{
				buf_ = new std::unique_ptr<double[]>[m_buffer_number];
				for(int i = 0; i < buffer_number; i++)
				{
					buf_[i]=std::unique_ptr<double[]>(new double[max_size_]);
				}
			};
		void Put(double const item[]);
		void Get(double returned_data[]);
		double Look(int index, int channel);
		void Dump(int marker, double* history, double* shot);
		void Reset();
		bool Empty() const;
		bool Full() const;
		size_t Capacity() const;
		size_t Size() const;
	private:
		std::mutex mutex_;
		std::unique_ptr<double[]> *buf_;
		size_t head_ = 0;
		size_t tail_ = 0;
		const size_t max_size_;
		bool full_ = 0;
		const size_t m_buffer_number;
};
#endif
