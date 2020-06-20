#include <math.h>

#include <gsl/gsl_blas.h>

#include "qdm.h"

#define TAU_INITIAL_GUESS 0.5
#define TAU_EPS 0.001
#define TAU_ITERATION_MAX 1000

#define RESET_SIZE 30
const double RESET_VALUES[RESET_SIZE] = {
  0.010000000000000000,
  0.010100000000000000,
  0.010200000000000001,
  0.010300000000000000,
  0.010400000000000000,
  0.010500000000000001,
  0.010600000000000000,
  0.010699999999999999,
  0.010800000000000001,
  0.010900000000000000,
  0.010999999999999999,
  0.010999999999999999,
  0.150714285714285717,
  0.290428571428571425,
  0.430142857142857160,
  0.569857142857142840,
  0.709571428571428520,
  0.849285714285714310,
  0.988999999999999990,
  0.988999999999999990,
  0.989099999999999979,
  0.989199999999999968,
  0.989299999999999957,
  0.989399999999999946,
  0.989500000000000046,
  0.989600000000000035,
  0.989700000000000024,
  0.989800000000000013,
  0.989900000000000002,
  0.989999999999999991,
};

/* Helper: Compute ispline by mmm dot product. */
static
int
ispline_mmm(
    double *result,

    double tau,
    size_t spline_df,
    const gsl_vector *knots,
    const gsl_vector *mmm
)
{
  int status = 0;

  size_t m = (knots->size - spline_df) + 1;
  gsl_vector *ispline = gsl_vector_alloc(m);

  qdm_ispline_vector(ispline, tau, spline_df, knots);
  status = gsl_blas_ddot(ispline, mmm, result);
  if (status != 0) {
    goto cleanup;
  }

cleanup:
  gsl_vector_free(ispline);

  return status;
}

/* Helper: Compute mspline by mmm dot product. */
static
int
mspline_mmm(
    double *result,

    double tau,
    size_t spline_df,
    const gsl_vector *knots,
    const gsl_vector *mmm
)
{
  int status = 0;

  size_t m = (knots->size - spline_df) + 1;
  gsl_vector *mspline = gsl_vector_alloc(m);

  qdm_mspline_vector(mspline, tau, spline_df, knots);
  status = gsl_blas_ddot(mspline, mmm, result);
  if (status != 0) {
    goto cleanup;
  }

cleanup:
  gsl_vector_free(mspline);

  return status;
}


int
qdm_find_tau(
    double *result,

    double v,
    size_t spline_df,
    const gsl_vector *knots,
    const gsl_vector *mmm
)
{
  int status = 0;

  double tau = TAU_INITIAL_GUESS;
  double qi_u = 0;
  double qm_u = 0;

  size_t reset_i = 0;
  for (size_t i = 0; i < TAU_ITERATION_MAX; i++) {
    if (reset_i >= RESET_SIZE) {
      tau = TAU_INITIAL_GUESS;

      break;
    }

    /* Calculate position... */
    status = ispline_mmm(&qi_u, tau, spline_df, knots, mmm);
    if (status != 0) {
      goto cleanup;
    }

    /* Calculate slope... */
    status = mspline_mmm(&qm_u, tau, spline_df, knots, mmm);
    if (status != 0) {
      goto cleanup;
    }

    /* Check if the update is within our desired interval [0, 1]. */
    double update = tau - (qi_u - v) / qm_u;
    if (update < 0 || update > 1) {
      tau = RESET_VALUES[reset_i];
      reset_i++;

      continue;
    } else {
      tau = update;
    }

    if (fabs(qi_u - v) <= TAU_EPS) {
      break;
    }
  }

  *result = tau;

cleanup:

  return status;
}

int
qdm_logl(
    double *log_likelihood,
    double *tau,

    double v,
    size_t spline_df,
    const gsl_vector *knots,
    const gsl_vector *mmm,

    double tau_low,
    double tau_high,
    double xi_low,
    double xi_high
)
{
  int status = 0;

  double low_threshold = 0;
  double high_threshold = 0;

  status = ispline_mmm(&low_threshold, tau_low, spline_df, knots, mmm);
  if (status != 0) {
    goto cleanup;
  }

  status = ispline_mmm(&high_threshold, tau_high, spline_df, knots, mmm);
  if (status != 0) {
    goto cleanup;
  }

  if (v <= low_threshold) {
    double sigma_low = 0;

    status = mspline_mmm(&sigma_low, tau_low, spline_df, knots, mmm);
    if (status != 0) {
      goto cleanup;
    }

    sigma_low *= tau_low;

    double z = low_threshold - v;

    *log_likelihood = log(tau_low / sigma_low * pow(1 + (xi_low * z) / sigma_low, -1 / xi_low - 1));
    *tau = (1 - (1 - pow(1 + xi_low * z / sigma_low, -1 / xi_low))) * tau_low;
  } else if (v >= high_threshold) {
    double sigma_high = 0;

    status = mspline_mmm(&sigma_high, tau_high, spline_df, knots, mmm);
    if (status != 0) {
      goto cleanup;
    }

    sigma_high *= (1 - tau_high);

    double z = v - high_threshold;

    *log_likelihood = log((1 - tau_high) / sigma_high * pow(1 + (xi_high * z) / sigma_high, -1 / xi_high - 1));
    *tau = (1 - pow(1 + xi_high * z / sigma_high, -1 / xi_high)) * (1 - tau_high) + tau_high;
  } else {
    double d = 0;

    status = qdm_find_tau(tau, v, spline_df, knots, mmm);
    if (status != 0) {
      goto cleanup;
    }

    status = mspline_mmm(&d, *tau, spline_df, knots, mmm);
    if (status != 0) {
      goto cleanup;
    }

    *log_likelihood = log(1 / d);
  }

cleanup:

  return status;
}

int
qdm_logl_2(
    double *log_likelihood,
    double *tau,

    double x,
    double y,

    double tau_low,
    double tau_high,

    double xi_low,
    double xi_high,

    size_t spline_df,

    const gsl_matrix *theta,
    const gsl_vector *knots
)
{
  int status = 0;

  double xi_data[2] = {1, x};
  gsl_vector_view xi = gsl_vector_view_array(xi_data, 2);

  double mmm_data[theta->size2];
  gsl_vector_view mmm = gsl_vector_view_array(mmm_data, theta->size2);

  status = gsl_blas_dgemv(CblasTrans, 1, theta, &xi.vector, 0, &mmm.vector);
  if (status != 0) {
    return status;
  }

  status = qdm_logl(
      log_likelihood,
      tau,

      y,
      spline_df,
      knots,
      &mmm.vector,

      tau_low,
      tau_high,
      xi_low,
      xi_high
  );
  if (status != 0) {
    return status;
  }

  return status;
}

void
qdm_logl_3(
    double *log_likelihood,
    double *tau,

    double x,
    double y,

    const qdm_tau *t,
    const gsl_vector *xi,

    const gsl_matrix *theta
)
{
  double xi_intercept_data[2] = {1, x};
  gsl_vector_view xi_intercept = gsl_vector_view_array(xi_intercept_data, 2);

  double mmm_data[theta->size2];
  gsl_vector_view mmm = gsl_vector_view_array(mmm_data, theta->size2);

  gsl_blas_dgemv(CblasTrans, 1, theta, &xi_intercept.vector, 0, &mmm.vector);

  double low_threshold = qdm_tau_ispline_mmm(t, t->low, &mmm.vector);
  double high_threshold = qdm_tau_ispline_mmm(t, t->high, &mmm.vector);

  /*
  fprintf(stderr, "x    : %f\n", x);
  fprintf(stderr, "y    : %f\n", y);

  fprintf(stderr, "tl   : %f\n", t->low);
  fprintf(stderr, "th   : %f\n", t->high);

  fprintf(stderr, "low  : %f\n", low_threshold);
  fprintf(stderr, "high : %f\n", high_threshold);

  qdm_matrix_csv_fwrite(stderr, theta);
  */

  double xi_low = gsl_vector_get(xi, 0);
  double xi_high = gsl_vector_get(xi, 1);

  if (y <= low_threshold) {
    double sigma_low = qdm_tau_mspline_mmm(t, t->low, &mmm.vector) * t->low;
    double z = low_threshold - y;

    *log_likelihood = log(t->low / sigma_low * pow(1 + (xi_low * z) / sigma_low, -1 / xi_low - 1));
    *tau = (1 - (1 - pow(1 + xi_low * z / sigma_low, -1 / xi_low))) * t->low;
  } else if (y >= high_threshold) {
    double sigma_high = qdm_tau_mspline_mmm(t, t->high, &mmm.vector) * (1 - t->high);
    double z = y - high_threshold;

    *log_likelihood = log((1 - t->high) / sigma_high * pow(1 + (xi_high * z) / sigma_high, -1 / xi_high - 1));
    *tau = (1 - pow(1 + xi_high * z / sigma_high, -1 / xi_high)) * (1 - t->high) + t->high;
  } else {
    *tau = qdm_tau_find(t, y, &mmm.vector);
    *log_likelihood = log(1 / qdm_tau_mspline_mmm(t, *tau, &mmm.vector));
  }

  /*
  fprintf(stderr, "ll   : %f\n", *log_likelihood);
  fprintf(stderr, "tau  : %f\n", *tau);
  */
}
