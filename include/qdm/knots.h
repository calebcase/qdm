#ifndef QDM_KNOTS_H
#define QDM_KNOTS_H 1

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

gsl_vector *
qdm_knots_vector(
  size_t spline_df,
  gsl_vector *knots_inter
);

int
qdm_knots_optimize(
  gsl_vector *result,

  gsl_rng *rng,

  gsl_vector *sorted_data,
  gsl_vector *middle,
  gsl_vector *possible_knots,

  size_t iterate_n,
  size_t spline_df
);

#endif /* QDM_KNOTS_H */
