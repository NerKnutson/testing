#ifndef NLMS_H
#define NLMS_H
class NLMS
{
	public:
	// Constructor (instantiates constants immediately) and then creates arrays
		NLMS(int N_chan, int filter_order, double step_size, double alpha, double tiny) :
			m_N_chan(N_chan),
			m_filter_order(filter_order),
			m_step_size(step_size),
			m_alpha(alpha),
			m_tiny(tiny),
			index(0)
			{
				m_x = new double*[N_chan];
				m_coeff = new double**[N_chan];
				m_error = new double[N_chan];
				m_var = new double[N_chan];
				m_triggered = new bool[N_chan];
				for (int i = 0; i < N_chan; i++)
				{
					m_triggered[i] = false;
					m_x[i] = new double[filter_order];
					m_coeff[i] = new double*[N_chan];
					for (int j = 0; j < N_chan; j++)
						m_coeff[i][j] = new double[filter_order + 1];
				}
			};
	// Destructor
		~NLMS()
			{
				delete m_x;
				delete m_coeff;
				delete m_error;
				delete m_var;
				delete m_triggered;
			};
		bool Process(double data[]);
		bool IsTriggered(int &which_channel);
		double* m_error;
		double* m_var;
		bool* m_triggered;
	private:
		const int m_N_chan;
		int index;
		int m_filter_order;
		double m_step_size;
		double m_alpha;
		const double m_tiny;
		double** m_x;
		double*** m_coeff;
};
#endif
