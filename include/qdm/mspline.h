#ifndef QDM_MSPLINE_H
#define QDM_MSPLINE_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void
qdm_mspline_vector(gsl_vector *result, const double tau, const size_t spline_df, const gsl_vector *knots);

void
qdm_mspline_matrix(gsl_matrix *result, const size_t spline_df, const gsl_vector *x, const gsl_vector *knots);

#endif /* QDM_MSPLINE_H */
