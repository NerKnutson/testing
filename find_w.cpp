#include <string.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "find_w.h"

// constants
//const gsl_complex zzero = {.dat[0] = 0.0, .dat[1] = 0.0};
//const gsl_complex zone = {.dat[0] = 1.0, .dat[1] = 0.0};
//const gsl_complex ztwo = {.dat[] = 2.0, .dat[1] = 0.0};
//const gsl_complex zmone = {.dat[0] = -1.0, .dat[1] = 0.0};
const gsl_complex zzero = {0.0, 0.0};
const gsl_complex zone = {1.0, 0.0};
const gsl_complex ztwo = {2.0, 0.0};
const gsl_complex zmone = {-1.0, 0.0};

// static function prototypes
static void calc_steering_vecs(const gsl_vector* slowness, opt_params* params);
static void calc_s_hat(opt_params* params);
static int calc_real_residual_vec(const gsl_vector* slowness, void* parameters,
				  gsl_vector* residual);
static int calc_real_jacobian(const gsl_vector* slowness, void* parameters,
			      gsl_matrix* jacobian);
static void wgf_fill_jacobian(opt_params* params, gsl_matrix* jacobian);
static void wgf_vector_complex_zero_the_imag_part(gsl_vector_complex* x);
static opt_params* opt_params_alloc(int n_fft, double sps, int n_omegas,
				    int* list_bins, gsl_matrix* R);

// functions
int init_find_w(const char* filename, opt_params** find_w_opt_params)
{

	// open the initialization file
  FILE* fp = fopen(filename, "r");

  // read in the parameters
  char str[10], array_filename[20];

  // FFT size
  int n_fft;
  if (fscanf(fp, "%s %d", str, &n_fft) != 2) return -1;
  if (strcmp("n_fft", str) != 0) return -1;

  // starting bin index for processing
  int i_start;
  if (fscanf(fp, "%s %d", str, &i_start) != 2) return -2;
  if (strcmp("i_start", str) != 0) return -2;

  // stopping bin index for processing
  int i_stop;
  if (fscanf(fp, "%s %d", str, &i_stop) != 2) return -3;
  if (strcmp("i_stop", str) != 0) return -3;

  // number of microphones
  int n_mics;
  if (fscanf(fp, "%s %d", str, &n_mics) != 2) return -4;
  if (strcmp("n_mics", str) != 0) return -4;

  // sampling rate
  double sps;
  if (fscanf(fp, "%s %lf", str, &sps) != 2) return -5;
  if (strcmp("sps", str) != 0) return -5;

  // parameter convergence tolerance
  double xtol;
  if (fscanf(fp, "%s %lf", str, &xtol) != 2) return -6;
  if (strcmp("xtol", str) != 0) return -6;

  // gradient tolerance
  double gtol;
  if (fscanf(fp, "%s %lf", str, &gtol) != 2) return -7;
  if (strcmp("gtol", str) != 0) return -7;

  // function convergence tolerance
  double ftol;
  if (fscanf(fp, "%s %lf", str, &ftol) != 2) return -8;
  if (strcmp("ftol", str) != 0) return -8;

  // max iterations
  int max_iter;
  if (fscanf(fp, "%s %d", str, &max_iter) != 2) return -9;
  if (strcmp("max_iter", str) != 0) return -9;

  // array geometry filename
  if (fscanf(fp, "%s", array_filename) != 1) return -10;

  // close file
  fclose(fp);

  // create the list of bins
  int n_omegas = i_stop - i_start + 1;
  int* list_bins = (int*) malloc(n_omegas*sizeof(int));
  int k = 0; for (int i = 0; i < n_omegas; ++i) list_bins[i] = k++;

  // read in the array geometry
  gsl_matrix* R = gsl_matrix_alloc(n_mics, 3);
  fp = fopen(array_filename, "r");
  double x;
  for (int i = 0; i < n_mics; ++i)
    for (int j = 0; j < 3; ++j) {
      if (fscanf(fp, "%lf", &x) != 1) return -11;
      gsl_matrix_set(R, i, j, x);
    }
  fclose(fp);

  // create an opt_params structure
  *find_w_opt_params = opt_params_alloc(n_fft, sps, n_omegas, list_bins, R);

  // assign remaining values
  (*find_w_opt_params)->maxiter = max_iter;
  (*find_w_opt_params)->xtol = xtol;
  (*find_w_opt_params)->gtol = gtol;
  (*find_w_opt_params)->ftol = ftol;

  // free memory
  gsl_matrix_free(R);
  free(list_bins);

  // return
  return 0;
}

/******************************************************************************/

void perform_opt(opt_params* params, gsl_matrix_complex** cps,
		 gsl_vector_complex** y, int mic1) {

  // index into bins to use
  int j;

  // calculate the Cholesky factorizations of the cross power spectra (cps),
  // l_inv_y, and the inverse of the lower factor. y are the fft data.
  for (int i = 0; i < params->n_omegas; ++i) {
    j = params->list_bins[i];

    // Cholesky
    gsl_matrix_complex_memcpy(params->l_inv[i],cps[j]);

    gsl_linalg_complex_cholesky_decomp(params->l_inv[i]);

    // invert lower triangle
    gsl_linalg_complex_tri_invert(CblasLower, CblasNonUnit, params->l_inv[i]);

    // get l_inv_y
    gsl_vector_complex_memcpy(params->l_inv_y[i], y[j]);
    gsl_blas_ztrmv(CblasLower, CblasNoTrans, CblasNonUnit, params->l_inv[i],
		   params->l_inv_y[i]);

    // now get L^-1 * minus_i_omega_R for use in jacobian calculations
    // first, make copy of minus i omega R
    //*params->l_inv_minus_i_omega_R[i] = *params->minus_i_omega_R[i];
    gsl_matrix_complex_memcpy(params->l_inv_minus_i_omega_R[i], params->minus_i_omega_R[i]);

    // now scale on the left by l_inv
    gsl_blas_ztrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
		   zone, params->l_inv[i], params->l_inv_minus_i_omega_R[i]);
  }

  // given the index of the mic that triggered this detection, calculate
  // an initial guess for the solution
  gsl_vector_view x0 = gsl_matrix_row(params->R, mic1);

  gsl_vector_scale(&x0.vector, -1.0/(343.0*gsl_blas_dnrm2(&x0.vector)));

// Uncomment to run ntk_test_routine(params)
//	calc_real_residual_vec(&x0.vector, params, params->residual);
//<@@>


  // initialize the optimizer for this run
  params->fdf.params = params;


  gsl_multifit_nlinear_init(&x0.vector, &params->fdf, params->w);

  // run till convergence
  params->status = gsl_multifit_nlinear_driver(params->maxiter, params->xtol,
					       params->gtol, params->ftol, NULL,
					       NULL, &params->info, params->w);

  // analyze results
  if (params->status == GSL_SUCCESS) {

    // residual vector and rss
    params->residual = gsl_multifit_nlinear_residual(params->w);
    double rss;
    gsl_blas_ddot(params->residual, params->residual, &rss);

    // unscaled covariance of the estimate
    gsl_matrix *J = gsl_multifit_nlinear_jac(params->w);
    params->covar = gsl_matrix_alloc (params->fdf.p, params->fdf.p);
    gsl_multifit_nlinear_covar(J, 1e-15, params->covar);

    // get scaled result
    double var = rss / (double)(params->fdf.n - params->fdf.p - 1);
    gsl_matrix_scale(params->covar, var);
  }
else if(params->status == GSL_ENOPROG)
{
	printf("GSL failed to find acceptable step\n");
	fflush(stdout);
}
else if(params->status == GSL_EMAXITER)
{
	printf("GSL failed to converge within maximimum iterations\n");
	fflush(stdout);
}
else
{
	printf("GSL straight up isn't doing it\n");
	fflush(stdout);
}

  // return
  // results are retained in the structure params
}

/******************************************************************************/

static void calc_steering_vecs(const gsl_vector* slowness, opt_params* params)
{
  gsl_vector_complex_view position;
  gsl_complex z;

  // make a complex slowness vector
  for (int i = 0; i < 3; ++i)
    GSL_SET_COMPLEX(gsl_vector_complex_ptr(params->zslowness, i),
		    gsl_vector_get(slowness, i), 0.0);

  // loop over frequencies to be used
  for (int i_omega = 0; i_omega < params->n_omegas; ++i_omega) {

    // loop over measurements [exp(-i * omega r \dot w)]
    for (int i_meas = 0; i_meas < params->n_meas; ++i_meas) {
      position = gsl_matrix_complex_row(params->minus_i_omega_R[i_omega],
					i_meas);
      gsl_blas_zdotu(&position.vector, params->zslowness, &z);
      gsl_vector_complex_set(params->steering_vecs[i_omega], i_meas,
			     gsl_complex_exp(z));
    }

    // now multiply by l_inv
    gsl_vector_complex_memcpy(params->l_inv_sv[i_omega], params->steering_vecs[i_omega]);
    gsl_blas_ztrmv(CblasLower, CblasNoTrans, CblasNonUnit,
		   params->l_inv[i_omega], params->l_inv_sv[i_omega]);

    // get the norm-squared
    gsl_blas_zdotc(params->l_inv_sv[i_omega], params->l_inv_sv[i_omega],
		   &params->norm_sq_l_inv_sv[i_omega]);
  }
}

/******************************************************************************/

static void calc_s_hat(opt_params* params)
{
  gsl_complex num;

  // loop over frequencies
  for (int i_omega = 0; i_omega < params->n_omegas; ++i_omega) {
    gsl_blas_zdotc(params->l_inv_sv[i_omega],
		   params->l_inv_y[i_omega],
		   &num);
    params->s_hat[i_omega] = gsl_complex_div(num,
			     params->norm_sq_l_inv_sv[i_omega]);
  }
}

/******************************************************************************/

static int calc_real_residual_vec(const gsl_vector* slowness, void* parameters,
				  gsl_vector* residual)
{
  opt_params* params = (opt_params*) parameters;

  // save slowness for reducting load on future jacobian calculations
  gsl_vector_memcpy(params->slowness,slowness);

  // calculate steering vectors
  calc_steering_vecs(slowness, params);

  // calculate s_hat
  calc_s_hat(params);

  int i_res = 0;
  gsl_complex z;

  // loop over frequencies
  for (int i_omega = 0; i_omega < params->n_omegas; ++i_omega) {

//    // loop over channels to get the weighted residual vector
#define Y(i, j) gsl_vector_complex_get(params->l_inv_y[(i)], (j))
#define Y_hat(i, j) gsl_complex_mul(params->s_hat[i], gsl_vector_complex_get(params->l_inv_sv[(i)], (j)))
    for (int i_meas = 0; i_meas < params->n_meas; ++i_meas) {
      z = gsl_complex_sub(Y(i_omega, i_meas), Y_hat(i_omega, i_meas));
      gsl_vector_set(residual, i_res++, GSL_REAL(z));
      gsl_vector_set(residual, i_res++, GSL_IMAG(z));
		//printf("%f + i %f\n",GSL_REAL(z), GSL_IMAG(z));
		//fflush(stdout);
    }
#undef Y
#undef Y_hat
  }

  return GSL_SUCCESS;
}

/******************************************************************************/

static int calc_real_jacobian(const gsl_vector* slowness, void* parameters,
			      gsl_matrix* jacobian)
{
  opt_params* params = (opt_params*) parameters;
  gsl_vector_complex_view row;
  gsl_complex z;

  // check if we can save some calculations
  if (gsl_vector_equal(slowness, params->slowness) == 0) {
    calc_steering_vecs(slowness, params);
    calc_s_hat(params);
  }

  // loop over frequencies
  for (int i_omega = 0; i_omega < params->n_omegas; ++i_omega) {

    // loop over channels to get partial of scaled sv wrt slowness
    for (int i_meas = 0; i_meas < params->n_meas; ++i_meas) {
      gsl_matrix_complex_memcpy(params->l_inv_p_sv_p_w[i_omega]
	, params->l_inv_minus_i_omega_R[i_omega]);
      row = gsl_matrix_complex_row(params->l_inv_p_sv_p_w[i_omega], i_meas);
      gsl_vector_complex_scale(&row.vector,
	        gsl_vector_complex_get(params->steering_vecs[i_omega], i_meas));
    }

    // now get partial of minus s_hat wrt to slowness
    gsl_blas_zgemv(CblasConjTrans, zone, params->l_inv_p_sv_p_w[i_omega],
		   params->l_inv_y[i_omega], zzero, params->work1);
    gsl_blas_zgemv(CblasConjTrans, zone, params->l_inv_p_sv_p_w[i_omega],
		   params->l_inv_sv[i_omega], zzero, params->work2);

    wgf_vector_complex_zero_the_imag_part(params->work2);
    gsl_vector_complex_scale(params->work2, ztwo);
    gsl_vector_complex_scale(params->work2, params->s_hat[i_omega]);
    gsl_vector_complex_sub(params->work1, params->work2);
    //z = gsl_complex_mul(zmone, params->norm_sq_l_inv_sv[i_omega]);
    //gsl_vector_complex_scale(params->p_s_hat_p_w,
			     //gsl_complex_inverse(z));
    gsl_vector_complex_scale(params->p_s_hat_p_w,
			     gsl_complex_inverse(params->norm_sq_l_inv_sv[i_omega]));

    // finally, lets get the jacobian matrix values using the indexing system
    // the rows as i_omega, i_meas, (real, imag)
    // -l_inv_p_sv_p_w * s_hat - l_inv_sv * p_s_hat_p_w
    gsl_matrix_complex_scale(params->l_inv_p_sv_p_w[i_omega],
			     params->s_hat[i_omega]);
    gsl_blas_zgeru(zone, params->l_inv_sv[i_omega], params->p_s_hat_p_w,
		   params->l_inv_p_sv_p_w[i_omega]);

// Noah's line
    gsl_matrix_complex_scale(params->l_inv_p_sv_p_w[i_omega],
			     zmone);
  }

  // a dreaded floating point copy. refactor one day
  wgf_fill_jacobian(params, jacobian);

  return GSL_SUCCESS;
}

/******************************************************************************/

static void wgf_vector_complex_zero_the_imag_part(gsl_vector_complex* x) {

/*
gsl_vector_view temp = gsl_vector_complex_imag(x);
gsl_vector_set_zero(&temp.vector);
  for (size_t i = 0; i < x->size; ++i)
	{
	printf("%f + i %f\n", GSL_REAL(gsl_vector_complex_get(x,i)),GSL_IMAG(gsl_vector_complex_get(x,i)));
	fflush(stdout);
	}
*/
  for (size_t i = 0; i < x->size; ++i)
    GSL_SET_IMAG(gsl_vector_complex_ptr(x, i), 0.0);
}

/******************************************************************************/

static void wgf_fill_jacobian(opt_params* params, gsl_matrix* jacobian) {

  int jac_row1 = 0;
  int jac_row2 = 1;
  int jac_col;
  gsl_complex* x;

  // some pointer arithmetic--yuk
  for (int i_omega = 0; i_omega < params->n_omegas; ++i_omega)
  {
    for (int i_meas = 0; i_meas < params->n_meas; ++i_meas) {
      jac_col = 0;
      x = gsl_matrix_complex_ptr(params->l_inv_p_sv_p_w[i_omega], i_meas, 0);
      gsl_matrix_set(jacobian, jac_row1, jac_col, x->dat[0]);
      gsl_matrix_set(jacobian, jac_row2, jac_col, x->dat[1]);
      x = gsl_matrix_complex_ptr(params->l_inv_p_sv_p_w[i_omega], i_meas, 1);
      gsl_matrix_set(jacobian, jac_row1, ++jac_col, x->dat[0]);
      gsl_matrix_set(jacobian, jac_row2, jac_col, x->dat[1]);
      x = gsl_matrix_complex_ptr(params->l_inv_p_sv_p_w[i_omega], i_meas, 2);
      gsl_matrix_set(jacobian, jac_row1, ++jac_col, x->dat[0]);
      gsl_matrix_set(jacobian, jac_row2, jac_col, x->dat[1]);
      jac_row1 += 2;
      jac_row2 += 2;
/*
      jac_row2 = jac_row1 + 2;
      jac_col = 0;
      x = gsl_matrix_complex_ptr(params->l_inv_p_sv_p_w[i_omega], i_meas, 0);
      gsl_matrix_set(jacobian, jac_row1, jac_col, x->dat[0]);
      gsl_matrix_set(jacobian, jac_row2, jac_col, x->dat[1]);
      x = gsl_matrix_complex_ptr(params->l_inv_p_sv_p_w[i_omega], i_meas, 1);
      gsl_matrix_set(jacobian, jac_row1, ++jac_col, x->dat[0]);
      gsl_matrix_set(jacobian, jac_row2, jac_col, x->dat[1]);
      x = gsl_matrix_complex_ptr(params->l_inv_p_sv_p_w[i_omega], i_meas, 2);
      gsl_matrix_set(jacobian, jac_row1, ++jac_col, x->dat[0]);
      gsl_matrix_set(jacobian, jac_row2, jac_col, x->dat[1]);
*/
    }
  }
}

/******************************************************************************/

opt_params* opt_params_alloc(int n_fft, double sps, int n_omegas,
			     int* list_bins, gsl_matrix* R)
// note each row of R is a position vector
{
  // allocate for output
  opt_params* params = (opt_params*) malloc(sizeof(opt_params));

  // frequency bin width
  double df = sps / (double) n_fft;

  // 2 * pi * df
  double twopidf = 2.0 * M_PI * df;

  // assign values
  params->n_omegas = n_omegas;
  params->n_meas = R->size1;
  params->R = gsl_matrix_alloc(R->size1, 3);
  gsl_matrix_memcpy(params->R, R);
  params->maxiter = 100;
  params->xtol = 1e-8;
  params->gtol = 1e-8;
  params->ftol = 0.0;
  params->n_fft = n_fft;
  params->list_bins = (int*) malloc(n_omegas*sizeof(int));
  params->omegas = (double*) malloc(n_omegas*sizeof(double));
  for (int i = 0; i < n_omegas; ++i) {
    params->list_bins[i] = list_bins[i];
    params->omegas[i] = twopidf * (double) list_bins[i];
  }

  // allocate structure members
  params->residual = gsl_vector_alloc(2 * n_omegas * params->n_meas);
  params->covar = gsl_matrix_alloc(3, 3);
  params->slowness = gsl_vector_alloc(3);
  params->zslowness = gsl_vector_complex_alloc(3);
  params->s_hat = (gsl_complex*) malloc(n_omegas*sizeof(gsl_complex));
  params->norm_sq_l_inv_sv = (gsl_complex*)
    malloc(n_omegas*sizeof(gsl_complex));
  params->work1 = gsl_vector_complex_alloc(3);
  params->p_s_hat_p_w = params->work1;
  params->work2 = gsl_vector_complex_alloc(3);
  params->minus_i_omega_R = (gsl_matrix_complex**)
    malloc(n_omegas*sizeof(gsl_matrix_complex*));
  params->l_inv_minus_i_omega_R = (gsl_matrix_complex**)
    malloc(n_omegas*sizeof(gsl_matrix_complex*));
  params->l_inv = (gsl_matrix_complex**)
    malloc(n_omegas*sizeof(gsl_matrix_complex*));
  params->l_inv_p_sv_p_w = (gsl_matrix_complex**)
    malloc(n_omegas*sizeof(gsl_matrix_complex*));
  params->l_inv_p_s_hat_p_w = (gsl_matrix_complex**)
    malloc(n_omegas*sizeof(gsl_matrix_complex*));
  params->steering_vecs = (gsl_vector_complex**)
    malloc(n_omegas*sizeof(gsl_vector_complex*));
  params->l_inv_y = (gsl_vector_complex**)
    malloc(n_omegas*sizeof(gsl_vector_complex*));
  params->l_inv_sv = (gsl_vector_complex**)
    malloc(n_omegas*sizeof(gsl_vector_complex*));
  for (int i = 0; i < n_omegas; ++i) {
    params->l_inv[i] = gsl_matrix_complex_alloc(params->n_meas,params->n_meas);
    params->l_inv_y[i] = gsl_vector_complex_alloc(params->n_meas);
    params->l_inv_sv[i] = gsl_vector_complex_alloc(params->n_meas);
    params->steering_vecs[i] = gsl_vector_complex_alloc(params->n_meas);
    params->minus_i_omega_R[i] = gsl_matrix_complex_calloc(params->n_meas, 3);
    params->l_inv_minus_i_omega_R[i] =
      gsl_matrix_complex_calloc(params->n_meas, 3);
    params->l_inv_p_sv_p_w[i] = gsl_matrix_complex_alloc(params->n_meas, 3);
    params->l_inv_p_s_hat_p_w[i] = gsl_matrix_complex_alloc(params->n_meas, 1);
  }

  // now, calculate minus_i_omega_R once for all n_omegas
  gsl_complex* z;
  for (int i_omega = 0; i_omega < n_omegas; ++i_omega)
    for (size_t i_row = 0; i_row < R->size1; ++i_row)
      for (size_t i_col = 0; i_col < R->size2; ++i_col) {
	z = gsl_matrix_complex_ptr(params->minus_i_omega_R[i_omega], i_row,
				   i_col);
	GSL_SET_COMPLEX(z, 0.0, -params->omegas[i_omega]
			* gsl_matrix_get(R, i_row, i_col));
      }

  // set up the solver type
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;

  // choose default parameters for the solver
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();

  // specify the objective function
  params->fdf.f = calc_real_residual_vec;
  params->fdf.df = calc_real_jacobian;
  params->fdf.fvv = NULL;
  params->fdf.n = 2 * params->n_omegas * params->n_meas;
  params->fdf.p = 3;

  // allocate the workspace for the solver
  params->w = gsl_multifit_nlinear_alloc(T, &fdf_params, params->fdf.n,
					 params->fdf.p);

  return params;
}


void ntk_test_routine(opt_params* params)
{
	gsl_vector_complex** test_y;
		test_y = new gsl_vector_complex*[params->n_fft/2 - 1];
		for (int i_omega = 0; i_omega < params->n_fft/2 - 1; ++i_omega)
			test_y[i_omega] = gsl_vector_complex_calloc(params->n_meas);

	gsl_matrix_complex** test_cps;
		test_cps = new gsl_matrix_complex*[params->n_fft/2 - 1];
// Sets noise to very low
		for (int i_omega = 0; i_omega < params->n_fft/2 - 1; ++i_omega)
		{
			test_cps[i_omega] = gsl_matrix_complex_calloc(params->n_meas,params->n_meas);
			for (int i_meas = 0; i_meas < params->n_meas; ++i_meas)
				gsl_matrix_complex_set(
					test_cps[i_omega],
					i_meas, i_meas,
					gsl_complex_rect(0.1,0.0)
				);
		}

	int mic1 = 0;
	gsl_vector_complex* slowness;
		slowness = gsl_vector_complex_calloc(3);
	gsl_vector_complex_set_zero(slowness);
	gsl_vector_complex_set(
		slowness,
		mic1,
		gsl_complex_rect(-1.0/343.0,0)
	);

	gsl_vector_view position;
	gsl_vector_complex* zposition;
		zposition = gsl_vector_complex_calloc(3);

	gsl_complex z;

  double df = 4800.0 / (double) params->n_fft;

  // 2 * pi * df
  double twopidf = 2.0 * M_PI * df;

// Sets the EXACT value for data
	for (int i_omega = 0; i_omega <  params->n_fft/2 - 1; ++i_omega) {
	 // loop over measurements [exp(-i * omega r \dot w)]
		for (int i_meas = 0; i_meas < params->n_meas; ++i_meas) {
			position = gsl_matrix_row(params->R, i_meas);
				for (int i = 0; i < 3; i++)
					gsl_vector_complex_set(
						zposition,
						i,
						{0.0,gsl_vector_get(&position.vector,i)}
					);
			gsl_blas_zdotu(zposition,slowness, &z);
			gsl_vector_complex_set(
				test_y[i_omega],
				i_meas,
				gsl_complex_exp(
					gsl_complex_mul_real(
						z,-(i_omega + 1)*twopidf
					)
				)
			);
		}
	}

	perform_opt(params, test_cps, test_y, mic1);
}
