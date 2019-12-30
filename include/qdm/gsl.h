#ifndef QDM_GSL_H
#define QDM_GSL_H 1

#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <osqp/osqp.h>

gsl_vector *
qdm_vector_seq(double from, double to, double by);

size_t
qdm_vector_search(const gsl_vector *v, double needle);

void
qdm_vector_csv_fwrite(FILE *f, gsl_vector *v);

gsl_vector *
qdm_vector_csv_fread(FILE *stream);

void
qdm_matrix_csv_fwrite(FILE *f, gsl_matrix *m);

int
qdm_matrix_tmm(gsl_matrix *m, gsl_matrix *result);

int
qdm_matrix_det_tmm(gsl_matrix *m, double *det);

gsl_vector *
qdm_vector_quantile(gsl_vector *data, gsl_vector *probs);

double
qdm_vector_rss(const gsl_vector *y, const gsl_vector *fx);

double
qdm_vector_sum(gsl_vector *v);

void
qdm_matrix_select_upper_triangle(gsl_matrix *m);

int
qdm_matrix_to_csc_matrix(csc **result, gsl_matrix *m);

#endif /* QDM_GSL_H */
