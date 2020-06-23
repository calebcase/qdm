#ifndef QDM_TAU_H
#define QDM_TAU_H 1

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

typedef struct {
  bool use_table;

  double low;
  double high;

  size_t spline_df;
  const gsl_vector *knots;

  gsl_vector *reset;
  gsl_matrix *ispline_table;
  gsl_matrix *mspline_table;
} qdm_tau;

qdm_tau *
qdm_tau_alloc(
    int use_table,

    double low,
    double high,

    size_t spline_df,
    const gsl_vector *knots
);

void
qdm_tau_free(
    qdm_tau *t
);

double
qdm_tau_ispline_mmm(
    const qdm_tau *t,
    const double value,
    const gsl_vector *mmm
);

double
qdm_tau_mspline_mmm(
    const qdm_tau *t,
    const double value,
    const gsl_vector *mmm
);

double
qdm_tau_find(
    const qdm_tau *t,
    const double v,
    const gsl_vector *mmm
);

#endif /* QDM_TAU_H */
