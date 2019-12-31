#ifndef QDM_THETA_H
#define QDM_THETA_H 1

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

int
qdm_theta_optimize(
    gsl_vector *result,
    gsl_vector *emperical_quantiles,
    gsl_matrix *ix
);

void
qdm_theta_matrix_constrain(
    gsl_matrix *theta,
    double min
);

#endif /* QDM_THETA_H */
