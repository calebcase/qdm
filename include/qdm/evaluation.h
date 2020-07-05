#ifndef QDM_EVALUATION_H
#define QDM_EVALUATION_H 1

#include "mcmc.h"
#include "parameters.h"

typedef struct {
  // Parameters
  qdm_parameters parameters;

  // RNG
  gsl_rng *rng;

  // Years (pre-scaling)
  double years_min;
  double years_max;

  // Value Bounding
  double lower_bound;
  double upper_bound;

  // Diagnostics
  double waic;
  double pwaic;
  double dic;
  double pd;

  // Results
  qdm_tau *t;

  qdm_mcmc *mcmc;

  gsl_matrix *theta_bar;
  gsl_matrix *theta_star_bar;

  gsl_matrix *xi_cov;
  gsl_matrix *theta_star_cov;

  double xi_low_bar;
  double xi_high_bar;

  gsl_matrix *theta_acc;
  gsl_vector *xi_acc;

  gsl_vector *m_knots;

  unsigned long int rng_seed;

  // Performance Metrics
  double elapsed;
} qdm_evaluation;

qdm_evaluation *
qdm_evaluation_new(const qdm_parameters *parameters);

void
qdm_evaluation_free(qdm_evaluation *e);

void
qdm_evaluation_fprint(
    FILE *f,
    const qdm_evaluation *e
);

int
qdm_evaluation_write(
    hid_t id,
    const qdm_evaluation *e
);

int
qdm_evaluation_read(
    hid_t id,
    qdm_evaluation **e
);

int
qdm_evaluation_run(
    qdm_evaluation *e,
    const gsl_matrix *future,
    size_t data_year_idx,
    size_t data_month_idx,
    size_t data_value_idx
);

#endif /* QDM_EVALUATION_H */
