#include <stdlib.h>

#include <gsl/gsl_statistics_double.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#include "qdm.h"

qdm_ijk *
qdm_ijk_alloc(
    size_t size1,
    size_t size2,
    size_t size3
)
{
  qdm_ijk *t = malloc(sizeof(qdm_ijk));

  t->size1 = size1;
  t->size2 = size2;
  t->size3 = size3;

  t->data = malloc(sizeof(double) * size1 * size2 * size3);
  t->owner = true;

  return t;
}

qdm_ijk *
qdm_ijk_calloc(
    size_t size1,
    size_t size2,
    size_t size3
)
{
  qdm_ijk *t = malloc(sizeof(qdm_ijk));

  t->size1 = size1;
  t->size2 = size2;
  t->size3 = size3;

  t->data = calloc(sizeof(double), size1 * size2 * size3);
  t->owner = true;

  return t;
}

void
qdm_ijk_free(
    qdm_ijk *t
)
{
  if (t == NULL) {
    return;
  }

  if (t->owner) {
    free(t->data);
  }

  free(t);
}

qdm_ijk_view
qdm_ijk_view_array(
    double *base,
    size_t size1,
    size_t size2,
    size_t size3
)
{
  qdm_ijk_view v = {
    .ijk = {
      .size1 = size1,
      .size2 = size2,
      .size3 = size3,

      .data = base,
      .owner = false,
    },
  };

  return v;
}

gsl_vector_view
qdm_ijk_get_k(
    qdm_ijk *t,
    size_t i,
    size_t j
)
{
  return gsl_vector_view_array_with_stride(
      &t->data[i * t->size2 + j],
      t->size1 * t->size2,
      t->size3
  );
}

gsl_matrix_view
qdm_ijk_get_ij(
    qdm_ijk *t,
    size_t k
)
{
  return gsl_matrix_view_array(
      &t->data[k * t->size1 * t->size2],
      t->size1,
      t->size2
  );
}

double
qdm_ijk_get(
    qdm_ijk *t,
    size_t i,
    size_t j,
    size_t k
)
{
  gsl_matrix_view ij = qdm_ijk_get_ij(t, k);

  return gsl_matrix_get(&ij.matrix, i, j);
}

gsl_matrix *
qdm_ijk_cov(
    qdm_ijk *t
)
{
  size_t d = t->size1 * t->size2;
  size_t dr = 0;

  gsl_matrix *r = gsl_matrix_calloc(d, d);

  for (size_t i0 = 0; i0 < t->size1; i0++) {
    for (size_t j0 = 0; j0 < t->size2; j0++) {
      gsl_vector_view k0 = qdm_ijk_get_k(t, i0, j0);

      /* TODO: When we have time (and if it matters), this is symmetric and we
       * don't need to do a bunch of these if we get a little clever about
       * copying the mirror'd results.
       */
      for (size_t ip = 0; ip < t->size1; ip++) {
        for (size_t jp = 0; jp < t->size2; jp++) {
          gsl_vector_view kp = qdm_ijk_get_k(t, ip, jp);

          double cov = gsl_stats_covariance(
              k0.vector.data,
              k0.vector.stride,
              kp.vector.data,
              kp.vector.stride,
              kp.vector.size
          );

          gsl_matrix_set(r, dr, ip * t->size2 + jp, cov);
        }
      }

      dr++;
    }
  }

  return r;
}

int
qdm_ijk_write(
    hid_t id,
    const char *name,
    const qdm_ijk *t
)
{
  int status = 0;

  hid_t datatype  = -1;
  hid_t dataspace = -1;
  hid_t dataset   = -1;

  if (t == NULL) {
    goto cleanup;
  }

  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);

  status = H5Tset_order(datatype, H5T_ORDER_LE);
  if (status != 0) {
    goto cleanup;
  }

  hsize_t dims[3] = {
    t->size3,
    t->size1,
    t->size2,
  };
  dataspace = H5Screate_simple(3, dims, NULL);
  if (dataspace < 0) {
    status = dataspace;
    goto cleanup;
  }

  dataset = H5Dcreate(id, name, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {
    status = dataset;
    goto cleanup;
  }

  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, t->data);
  if (status != 0) {
    goto cleanup;
  }

cleanup: 
  if (dataset >= 0) {
    H5Dclose(dataset);
  }

  if (dataspace >= 0) {
    H5Sclose(dataspace);
  }

  if (datatype >= 0) {
    H5Tclose(datatype);
  }

  H5Oflush(id);

  return status;
}

int
qdm_ijk_read(
    hid_t id,
    const char *name,
    qdm_ijk **t
)
{
  int status = 0;

  int rank = 0;

  status = H5LTget_dataset_ndims(id, name, &rank);
  if (status < 0) {
    return status;
  }

  if (rank != 3) {
    return -1;
  }

  hsize_t dims[3];

  status = H5LTget_dataset_info(id, name, dims, NULL, NULL);
  if (status < 0) {
    return status;
  }

  size_t size3 = dims[0];
  size_t size1 = dims[1];
  size_t size2 = dims[2];

  qdm_ijk *tmp = qdm_ijk_alloc(size1, size2, size3);

  status = H5LTread_dataset_double(id, name, tmp->data);
  if (status < 0) {
    qdm_ijk_free(tmp);

    return status;
  }

  *t = tmp;

  return status;
}
