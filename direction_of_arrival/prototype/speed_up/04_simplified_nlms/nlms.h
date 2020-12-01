#ifndef NLMS_H
#define NLMS_H
class NLMS
{
	public:
	// Constructor (instantiates constants immediately) and then creates arrays
		NLMS(int N_chan, int filter_order, double step_size, double alpha, double threshold, double tiny, int head_channel) :
			m_N_chan(N_chan),
			m_filter_order(filter_order),
			m_step_size(step_size),
			m_alpha(alpha),
			m_threshold(threshold),
			m_tiny(tiny),
			magnitude(tiny),
			index(0),
				m_error(0),
				m_var(0),
				m_triggered(false),
				m_head_channel(head_channel)
			{
				m_x = new double*[N_chan];
				m_coeff = new double*[N_chan];
				for (int i = 0; i < N_chan; i++)
				{
					m_x[i] = new double[filter_order];
					m_coeff[i] = new double[filter_order + 1];
				}
			};
	// Destructor
		~NLMS()
			{
				for (int i = 0; i < m_N_chan; i++)
				{
					delete m_x[i];
					delete m_coeff[i];
				}
				delete m_x;
				delete m_coeff;
			};
	// Takes in a data array from N_chan channels
		bool Process(double data[]);
	// Returns true if triggered and requires an int to be passed to send back which channel was triggered
		bool IsTriggered(int &which_channel);
	// Returns true if triggered and requires an int to be passed to send back which channel was triggered
		void ResetTriggers();
	// Public member variables for each channel
		double m_error;
		double m_var;
		bool m_triggered;
	private:
	// Number of channels
		const int m_N_chan;
	// Magnitude of data vector
		double magnitude;
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
		double** m_coeff;
	// Channel that detects
		const int m_head_channel;
};
#endif
