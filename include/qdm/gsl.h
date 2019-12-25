#ifndef QDM_GSL_H
#define QDM_GSL_H 1

#include <stdio.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

gsl_vector *
qdm_vector_seq(double from, double to, double by);

size_t
qdm_vector_search(const gsl_vector *v, double needle);

void
qdm_vector_csv_fwrite(FILE *f, gsl_vector *v);

void
qdm_matrix_csv_fwrite(FILE *f, gsl_matrix *m);

#endif /* QDM_GSL_H */
