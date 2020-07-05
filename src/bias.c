#include <gsl/gsl_statistics_double.h>

#include "qdm.h"

gsl_matrix *
qdm_bias_values(
    qdm_evaluation *h,
    qdm_evaluation *f
)
{
  size_t n = f->mcmc->r.tau->size2;
  size_t s = f->mcmc->r.s;

  double tau_low = h->parameters.tau_low;
  double tau_high = h->parameters.tau_high;

  double lower_bound = h->lower_bound;
  double upper_bound = h->upper_bound;

  gsl_matrix *r = gsl_matrix_calloc(n, s);

  fprintf(stderr, "r[%zu, %zu]\n",
      r->size1,
      r->size2
  );

  fprintf(stderr, "r.tau[%zu, %zu, %zu]\n",
      f->mcmc->r.tau->size1,
      f->mcmc->r.tau->size2,
      f->mcmc->r.tau->size3
  );

  fprintf(stderr, "r.theta[%zu, %zu, %zu]\n",
      f->mcmc->r.theta->size1,
      f->mcmc->r.theta->size2,
      f->mcmc->r.theta->size3
  );

  fprintf(stderr, "r.xi[%zu, %zu, %zu]\n",
      f->mcmc->r.xi->size1,
      f->mcmc->r.xi->size2,
      f->mcmc->r.xi->size3
  );

  for (size_t j = 0; j < n; j++) {
    for (size_t k = 0; k < s; k++) {
      double tau = qdm_ijk_get(f->mcmc->r.tau, 0, j, k);
      double value = 0;

      gsl_matrix_view theta = qdm_ijk_get_ij(h->mcmc->r.theta, k);
      gsl_vector_view mmm = gsl_matrix_row(&theta.matrix, 0);

      if (tau_low <= tau && tau <= tau_high) {
        value = qdm_tau_ispline_mmm(
            h->t,
            tau,
            &mmm.vector
        );
      } else if (tau < tau_low) {
        double is = qdm_tau_ispline_mmm(
            h->t,
            tau_low,
            &mmm.vector
        );

        double ms = qdm_tau_ispline_mmm(
            h->t,
            tau_low,
            &mmm.vector
        );

        double xi_low = qdm_ijk_get(h->mcmc->r.xi, 0, 0, k);

        value = is - tau_low / ms * (pow(tau / tau_low, -xi_low) - 1);
        value = fmax(value, lower_bound);
      } else { // tau_high < tau
        double is = qdm_tau_ispline_mmm(
            h->t,
            tau_high,
            &mmm.vector
        );

        double ms = qdm_tau_ispline_mmm(
            h->t,
            tau_high,
            &mmm.vector
        );

        double xi_high = qdm_ijk_get(h->mcmc->r.xi, 0, 1, k);

        value = is - tau_high / ms * (pow((1 - tau) / (1 - tau_high), -xi_high) - 1);
        value = fmin(value, upper_bound);
      }

      gsl_matrix_set(r, j, k, value);
    }
  }

  return r;
}

gsl_matrix *
qdm_bias_correct(
    const gsl_vector *years,
    const double month,
    const gsl_vector *days,

    const gsl_vector *y,

    qdm_evaluation *fp,
    qdm_evaluation *ho,
    qdm_evaluation *hc
)
{
  gsl_matrix *hov = qdm_bias_values(ho, fp);
  gsl_matrix *hcv = qdm_bias_values(hc, fp);
  gsl_vector *row = gsl_vector_calloc(hov->size2);

  gsl_matrix *r = gsl_matrix_calloc(hov->size1, 5);

  fprintf(stderr, "years[%zu]\n", years->size);
  fprintf(stderr, "month %f\n", month);
  fprintf(stderr, "days[%zu]\n", days->size);
  fprintf(stderr, "y[%zu]\n", y->size);
  fprintf(stderr, "hov[%zu, %zu]\n", hov->size1, hov->size2);
  fprintf(stderr, "hcv[%zu, %zu]\n", hcv->size1, hcv->size2);
  fprintf(stderr, "row[%zu]\n", row->size);

  for (size_t i = 0; i < hov->size1; i++) {
    for (size_t j = 0; j < hov->size2; j++) {
      gsl_vector_set(row, j, gsl_matrix_get(hov, i, j) + (gsl_vector_get(y, i) - gsl_matrix_get(hcv, i, j)));

      double mean = gsl_stats_mean(row->data, row->stride, row->size);
      double sd = gsl_stats_sd(row->data, row->stride, row->size);

      gsl_matrix_set(r, i, 0, gsl_vector_get(years, i));
      gsl_matrix_set(r, i, 1, month);
      gsl_matrix_set(r, i, 2, gsl_vector_get(days, i));
      gsl_matrix_set(r, i, 3, mean);
      gsl_matrix_set(r, i, 4, sd);
    }
  }

  gsl_vector_free(row);
  gsl_matrix_free(hcv);
  gsl_matrix_free(hov);

  return r;
}
