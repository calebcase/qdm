#ifndef QDM_LOGL_H
#define QDM_LOGL_H 1

#include <gsl/gsl_vector.h>

int
qdm_find_tau(
    double *result,

    double v,
    size_t spline_df,
    const gsl_vector *knots,
    const gsl_vector *mmm
);

int
qdm_logl(
    double *log_likelihood,
    double *tau,

    double v,
    size_t spline_df,
    const gsl_vector *knots,
    const gsl_vector *mmm,

    double tau_low,
    double tau_high,
    double xi_low,
    double xi_high
);

#endif /* QDM_LOGL_H */
