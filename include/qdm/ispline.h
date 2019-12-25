#ifndef QDM_ISPLINE_H
#define QDM_ISPLINE_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void
qdm_ispline_vector(gsl_vector *result, const double tau, const size_t spline_df, const gsl_vector *knots);

void
qdm_ispline_matrix(gsl_matrix *result, const size_t spline_df, const gsl_vector *x, const gsl_vector *knots);

#endif /* QDM_ISPLINE_H */
