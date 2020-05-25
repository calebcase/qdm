#ifndef QDM_EVALUATION_H
#define QDM_EVALUATION_H 1

#include "mcmc.h"
#include "parameters.h"

typedef struct {
  // Parameters
  qdm_parameters parameters;

  // Years (pre-scaling)
  double years_min;
  double years_max;

  // Diagnostics
  double waic;
  double pwaic;
  double dic;
  double pd;

  // Results
  qdm_mcmc *mcmc;

  gsl_matrix *theta_bar;
  double xi_low_bar;
  double xi_high_bar;

  gsl_matrix *theta_acc;
  gsl_vector *xi_acc;

  gsl_vector *m_knots;

  gsl_matrix *bias_corrected;

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
qdm_evaluation_run(
    qdm_evaluation *e,
    const gsl_matrix *future,
    size_t data_year_idx,
    size_t data_month_idx,
    size_t data_value_idx
);

#endif /* QDM_EVALUATION_H */
