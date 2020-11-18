#ifndef RBUFFER_H
#define RBUFFER_H
#include <mutex>
#include <memory>
#include <cstddef>

class RingBuffer
{
	public:
	// Constructor sets the constants immediately and then takes care of the buffer declaration: buf[N_chan][buffer_size]
		explicit RingBuffer(size_t buffer_number, size_t size, int filter_order, double step_size, double alpha, double threshold, double tiny):
			m_buffer_number(buffer_number),
			max_size_(size),

			m_filter_order(filter_order),
			m_step_size(step_size),
			m_alpha(alpha),
			m_threshold(threshold),
			m_tiny(tiny),
			index(0)
			{
				buf_ = new std::unique_ptr<double[]>[m_buffer_number];
				for(int i = 0; i < m_buffer_number; i++)
				{
					buf_[i]=std::unique_ptr<double[]>(new double[max_size_]);
				}
				m_x = new double*[buffer_number];
				m_coeff = new double**[buffer_number];
				m_error = new double[buffer_number];
				m_var = new double[buffer_number];
				m_triggered = new bool[buffer_number];

				m_record_history = new bool[buffer_number];
				m_var_recorded = new bool[buffer_number];
				m_triggered_var = new double[buffer_number];

				for (int i = 0; i < buffer_number; i++)
				{
					m_triggered[i] = false;
					m_record_history[i] = true;
					m_var_recorded[i] = false;

					m_x[i] = new double[filter_order];
					m_coeff[i] = new double*[buffer_number];
					for (int j = 0; j < buffer_number; j++)
						m_coeff[i][j] = new double[filter_order + 1];
				}
			};
		~RingBuffer()
			{
				for(int i = 0; i < m_buffer_number; i++)
					buf_[i].reset();
				delete[] buf_;


				for (int i = 0; i < m_buffer_number; i++)
				{
					delete m_x[i];
					for (int j = 0; j < m_buffer_number; j++)
						delete m_coeff[i][j];
				}
				delete m_x;
				delete m_coeff;
				delete m_error;
				delete m_var;
				delete m_triggered;
				delete m_record_history;
				delete m_var_recorded;
				delete m_triggered_var;
			}
		void Put(double const item[]);
		void UnPut(int size);
		void Get(double returned_data[]);
		double Look(int index, int channel);
		void Dump(int marker, double* history, double* shot);
		void Reset();
		bool Empty() const;
		bool Full() const;
		size_t Capacity() const;
		size_t Size() const;


	// Takes in a data array from N_chan channels
		bool Process(double data[]);
	// Returns true if triggered and requires an int to be passed to send back which channel was triggered
		bool IsTriggered(int &which_channel);
	// Returns true if triggered and requires an int to be passed to send back which channel was triggered
		void ResetTriggers();
		bool IsRecording();
	// Public member variables for each channel
		double* m_error;
		double* m_var;
		bool* m_triggered;
		bool* m_record_history;
		bool* m_var_recorded;
		double* m_triggered_var;
	private:
		std::mutex mutex_;
		std::unique_ptr<double[]> *buf_;
		size_t head_ = 0;
		size_t tail_ = 0;
		const size_t max_size_;
		bool full_ = 0;
		const size_t m_buffer_number;


	// Dummy index for when using the Process function
		int index;
	// Number of data points in the past required from each channel
		int m_filter_order;
	// How quickly the algorithm is allowed to change AKA learning rate
		double m_step_size;
	// How much the past variance influences the current variance (usually close to 1)
		double m_alpha;
	// Sets threshold for detection
		double m_threshold;
	// Small number to not divide by zero. Possibly hardcode this in later?
		const double m_tiny;
	// Storing data
		double** m_x;
	// Coefficients for the algorithm
		double*** m_coeff;
};
#endif
