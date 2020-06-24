#ifndef QDM_IJK_H
#define QDM_IJK_H 1

#include <stddef.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <hdf5.h>

typedef struct {
  size_t size1;
  size_t size2;
  size_t size3;

  double *data;
  bool owner;
} qdm_ijk;

typedef struct {
  qdm_ijk ijk;
} qdm_ijk_view;

qdm_ijk *
qdm_ijk_alloc(
    size_t size1,
    size_t size2,
    size_t size3
);

qdm_ijk *
qdm_ijk_calloc(
    size_t size1,
    size_t size2,
    size_t size3
);

void
qdm_ijk_free(qdm_ijk *t);

qdm_ijk_view
qdm_ijk_view_array(
    double *base,
    size_t size1,
    size_t size2,
    size_t size3
);

gsl_vector_view
qdm_ijk_get_k(
    qdm_ijk *t,
    size_t i,
    size_t j
);

gsl_matrix_view
qdm_ijk_get_ij(
    qdm_ijk *t,
    size_t k
);

double
qdm_ijk_get(
    qdm_ijk *t,
    size_t i,
    size_t j,
    size_t k
);

gsl_matrix *
qdm_ijk_cov(
    qdm_ijk *t
);

int
qdm_ijk_write(
    hid_t id,
    const char *name,
    const qdm_ijk *t
);

#endif /* QDM_IJK_H */
