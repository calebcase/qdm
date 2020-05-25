#ifndef QDM_PARAMETERS_H
#define QDM_PARAMETERS_H 1

#include <gsl/gsl_vector.h>

typedef struct {
  int rng_seed;

  int acc_check;
  int burn;
  int iter;
  int knot_try;
  int spline_df;
  int thin;

  double month;
  double knot;

  double tau_high;
  double tau_low;

  double theta_min;
  double theta_tune_sd;

  double xi_high;
  double xi_low;
  double xi_prior_mean;
  double xi_prior_var;
  double xi_tune_sd;

  bool truncate;
} qdm_parameters;

void
qdm_parameters_fprint(
    FILE *f,
    const qdm_parameters *p
);

int
qdm_parameters_write(
    hid_t id,
    const qdm_parameters *p
);

#endif /* QDM_PARAMETERS_H */
