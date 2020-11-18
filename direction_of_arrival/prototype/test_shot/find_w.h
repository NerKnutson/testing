#ifndef FIND_W_H
#define FIND_W_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>

typedef struct
{
  int n_meas; // number of sensors
  int n_fft; // number of frequencies in fft
  int n_omegas; // number of frequences to use
  int* list_bins; // list of bins to use
  double* omegas; // radian frequencies used

  int status; // solver result
  int info; // solver result
  int maxiter; // max iterations
  double xtol; // delta x tolerance
  double gtol; // gradient tolerance;
  double ftol; // function value tolerance

  gsl_vector* slowness; // slowness vector
  gsl_vector_complex* zslowness; // slowness vector
  gsl_vector* residual; // residual vector
  gsl_matrix* covar; // slowness covariance estimate
  gsl_matrix* R; // matrix of mic position vectors (row wise)
  gsl_complex* s_hat; // signal fourier estimates
  gsl_complex* norm_sq_l_inv_sv; // norm-squared of l_inv_sv
  gsl_vector_complex* p_s_hat_p_w; // this is a pointer to work1
  gsl_vector_complex* work1; // a work vector
  gsl_vector_complex* work2; // a work vector
  gsl_vector_complex** steering_vecs; // steering vectors
  gsl_vector_complex** l_inv_y; // inverse of cholesky factors times ffts
  gsl_vector_complex** l_inv_sv; // inverse of cholesky factors times steering
  gsl_vector_complex** y; // a pointer to the measurements
  gsl_matrix_complex** minus_i_omega_R; // omegas times sensor position vectors
  gsl_matrix_complex** l_inv_minus_i_omega_R; // scaled above
  gsl_matrix_complex** l_inv; // inverses of cholesky factors of cross pwr spec
  gsl_matrix_complex** l_inv_p_sv_p_w; // partial wrt slowness
  gsl_matrix_complex** l_inv_p_s_hat_p_w; // another one

  gsl_multifit_nlinear_workspace *w; // optimization workspace pointer
  gsl_multifit_nlinear_fdf fdf; // optimization interface to the cost function

} opt_params;

int init_find_w(const char* filename, opt_params** find_w_opt_params);

void perform_opt(opt_params* params, gsl_matrix_complex** cps,
		 gsl_vector_complex** y, int mic1);

void ntk_test_routine(opt_params* params);

#endif
