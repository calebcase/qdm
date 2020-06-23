#ifndef QDM_BIAS_H
#define QDM_BIAS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "evaluation.h"

gsl_matrix *
qdm_bias_values(
    qdm_evaluation *h,
    qdm_evaluation *f
);

gsl_matrix *
qdm_bias_correct(
    const gsl_vector *years,
    const double month,
    const gsl_vector *days,

    const gsl_vector *y,

    qdm_evaluation *fp,
    qdm_evaluation *ho,
    qdm_evaluation *hc
);

#endif /* QDM_BIAS_H */
