#ifndef QDM_GSL_H
#define QDM_GSL_H 1

#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <osqp/osqp.h>

gsl_vector *
qdm_vector_seq(double from, double to, double by);

void
qdm_vector_set_seq(gsl_vector *v, double from, double to);

size_t
qdm_vector_search(const gsl_vector *v, double needle);

void
qdm_vector_csv_fwrite(FILE *f, const gsl_vector *v);

void
qdm_matrix_csv_fwrite(FILE *f, const gsl_matrix *m);

int
qdm_matrix_tmm(gsl_matrix *m, gsl_matrix *result);

int
qdm_matrix_det_tmm(gsl_matrix *m, double *det);

gsl_vector *
qdm_vector_sorted(const gsl_vector *v);

gsl_vector *
qdm_vector_quantile(gsl_vector *data, gsl_vector *probs);

double
qdm_vector_rss(const gsl_vector *y, const gsl_vector *fx);

double
qdm_vector_sum(gsl_vector *v);

size_t
qdm_vector_greater_than(const gsl_vector *v, double value);

void
qdm_matrix_select_upper_triangle(gsl_matrix *m);

int
qdm_matrix_to_csc_matrix(csc **result, gsl_matrix *m);

gsl_vector *
qdm_vector_copy(const gsl_vector *src);

gsl_matrix *
qdm_matrix_copy(const gsl_matrix *src);

int
qdm_vector_hd5_read(
    hid_t id,
    const char *name,
    gsl_vector **v
);

int
qdm_vector_hd5_write(
    hid_t id,
    const char *name,
    const gsl_vector *v
);

int
qdm_matrix_hd5_read(
    hid_t id,
    const char *name,
    gsl_matrix **m
);

int
qdm_matrix_hd5_write(
    hid_t id,
    const char *name,
    const gsl_matrix *m
);

gsl_vector *
qdm_matrix_filter(
    const gsl_matrix *m,
    size_t needle_column,
    double needle_value,
    size_t select_column
);

#endif /* QDM_GSL_H */
