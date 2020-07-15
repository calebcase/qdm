#ifndef QDM_CONFIG_H
#define QDM_CONFIG_H 1

#include <gsl/gsl_vector.h>

typedef struct {
  unsigned long int rng_seed;

  int acc_check;

  int burn_discovery;
  int iter_discovery;

  int burn_analysis;
  int iter_analysis;

  int knot_try;
  int spline_df;
  int thin;

  gsl_vector *months;
  gsl_vector *knots;

  gsl_vector *tau_highs;
  gsl_vector *tau_lows;

  double theta_min;
  double theta_tune_sd;

  double xi_high;
  double xi_low;
  double xi_prior_mean;
  double xi_prior_var;
  double xi_tune_sd;

  double bound;

  int debug;
  int tau_table;
} qdm_config;

int
qdm_config_read(
    qdm_config **config,
    const char *file_path,
    const char *group_path
);

void
qdm_config_fwrite(
    FILE *f,
    const qdm_config *cfg
);

#endif /* QDM_CONFIG_H */
