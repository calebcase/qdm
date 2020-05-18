#ifndef QDM_LOGL_H
#define QDM_LOGL_H 1

#include <gsl/gsl_matrix.h>
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

int
qdm_logl_2(
    double *log_likelihood,
    double *tau,

    double x,
    double y,

    double tau_low,
    double tau_high,

    double xi_low,
    double xi_high,

    size_t spline_df,

    const gsl_matrix *theta,
    const gsl_vector *knots
);

#endif /* QDM_LOGL_H */
