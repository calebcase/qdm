#ifndef QDM_KNOTS_H
#define QDM_KNOTS_H 1

#include <gsl/gsl_vector.h>

gsl_vector *
qdm_knots_vector(size_t spline_df, gsl_vector *knots_inter);

#endif /* QDM_KNOTS_H */
