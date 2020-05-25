/* Utility functions for working with GSL data types. */

#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_statistics.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <osqp/osqp.h>

#include "qdm.h"

gsl_vector *
qdm_vector_seq(double from, double to, double by)
{
  size_t size = fabs(to - from) / by;
  gsl_vector *s = gsl_vector_alloc(size + 1);

  double value = from;
  for (size_t i = 0; i < s->size; i++) {
    gsl_vector_set(s, i, value);
    value += by;
  }

  return s;
}

void
qdm_vector_set_seq(gsl_vector *v, double from, double to)
{
  double by = (to - from) / (double)(v->size - 1);

  double value = from;
  for (size_t i = 0; i < v->size; i++) {
    gsl_vector_set(v, i, value);
    value += by;
  }
}

size_t
qdm_vector_search(const gsl_vector *v, double needle)
{
  double x = 0;
  size_t i = v->size - 1;

  for (; i > 0; i--) {
    x = gsl_vector_get(v, i);
    
    if (needle >= x) {
      break;
    }
  }

  return i;
}

void
qdm_vector_csv_fwrite(FILE *f, const gsl_vector *v)
{
  for (size_t i = 0; i < v->size; i++) {
    fprintf(f, "%.17g", gsl_vector_get(v, i));
    if (i < v->size - 1) {
      fprintf(f, ",");
    }
  }
  fprintf(f, "\n");
}

gsl_vector *
qdm_vector_csv_fread(FILE *stream)
{
  const size_t chunk_size = 1024;
  size_t chunk_count = 0;

  size_t data_count = 0;
  double *data = NULL;

  int nread = 0;
  double value = 0;
  while (1) {
    nread = fscanf(stream, "%lg\n", &value);
    if (nread <= 0) {
      break;
    }

    if (data_count >= chunk_count * chunk_size) {
      chunk_count++;
      data = realloc(data, chunk_count * chunk_size * sizeof(double));
    }

    data[data_count] = value;
    data_count++;
  }

  gsl_vector *v = gsl_vector_alloc(data_count);
  for (size_t i = 0; i < data_count; i++) {
    gsl_vector_set(v, i, data[i]);
  }

  free(data);

  return v;
}

void
qdm_matrix_csv_fwrite(FILE *f, const gsl_matrix *m)
{
  if (m != NULL) {
    for (size_t i = 0; i < m->size1; i++) {
      gsl_vector_const_view row = gsl_matrix_const_row(m, i);
      qdm_vector_csv_fwrite(f, &row.vector);
    }
  }
}

/* Compute M^T * M */
int
qdm_matrix_tmm(gsl_matrix *m, gsl_matrix *result)
{
  return gsl_blas_dgemm(
      CblasTrans , CblasNoTrans , 1.0 ,
      m          , m            ,
      0.0        , result
  );
}

/* Compute det(M^T * M) */
int
qdm_matrix_det_tmm(gsl_matrix *m, double *det)
{
  int status = 0;
  gsl_matrix *c = gsl_matrix_alloc(m->size2, m->size2);
  gsl_permutation *p = gsl_permutation_alloc(c->size1);

  status = gsl_blas_dgemm(
      CblasTrans , CblasNoTrans , 1.0 ,
      m          , m            ,
      0.0        , c
  );
  if (status != 0) {
    goto cleanup;
  }

  int signum = 0;
  status = gsl_linalg_LU_decomp(c, p, &signum);
  if (status != 0) {
    goto cleanup;
  }

  *det = gsl_linalg_LU_det(c, signum);

cleanup:
  gsl_permutation_free(p);
  gsl_matrix_free(c);

  return status;
}

/* Create a sorted copy of the vector v. */
gsl_vector *
qdm_vector_sorted(const gsl_vector *v)
{
  gsl_vector *s = gsl_vector_alloc(v->size);

  gsl_vector_memcpy(s, v);
  gsl_sort_vector(s);

  return s;
}

gsl_vector *
qdm_vector_quantile(gsl_vector *data, gsl_vector *probs)
{
  gsl_vector *sorted = qdm_vector_sorted(data);
  gsl_vector *quantiles = gsl_vector_alloc(probs->size);

  for (size_t i = 0; i < probs->size; i++) {
    double q = gsl_stats_quantile_from_sorted_data(
        sorted->data,
        sorted->stride,
        sorted->size,
        gsl_vector_get(probs, i)
    );

    gsl_vector_set(quantiles, i, q);
  }

  free(sorted);

  return quantiles;
}

/* Compute the residual sum of squares:
 *
 * sum((y[i] - f(x[i])) ^ 2)
 */
double
qdm_vector_rss(const gsl_vector *y, const gsl_vector *fx)
{
  int status = 0;
  double rss = 0.0;

  gsl_vector *se = gsl_vector_alloc(y->size);
  status = gsl_vector_memcpy(se, y);
  if (status != 0) {
    goto cleanup;
  }

  status = gsl_vector_sub(se, fx);
  if (status != 0) {
    goto cleanup;
  }

  status = gsl_vector_mul(se, se);
  if (status != 0) {
    goto cleanup;
  }

  rss = qdm_vector_sum(se);

cleanup:
  gsl_vector_free(se);

  return rss;
}

/* Compute the summation of the vector.
 *
 * This uses the "iterative Kahan-Babuska algorithm" (aka Klein summation) as
 * outlined here:
 *
 * https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 */
double
qdm_vector_sum(gsl_vector *v)
{
  double s = 0.0;
  double cs = 0.0;
  double ccs = 0.0;

  for (size_t i = 0; i < v->size; i++) {
    double t = s + gsl_vector_get(v, i);

    double c;
    if (fabs(s) >= fabs(gsl_vector_get(v, i))) {
      c = (s - t) + gsl_vector_get(v, i);
    } else {
      c = (gsl_vector_get(v, i) - t) + s;
    }

    s = t;
    t = cs + c;

    double cc;
    if (fabs(cs) >= fabs(c)) {
      cc = (cs - t) + c;
    } else {
      cc = (c - t) + cs;
    }

    cs = t;
    ccs = ccs + cc;
  }

  return s + cs + ccs;
}

/* Return the first index with a value greater than the value. */
size_t
qdm_vector_greater_than(const gsl_vector *v, double value)
{
  for (size_t i = 0; i < v->size; i++) {
    if (gsl_vector_get(v, i) > value) {
      return i;
    }
  }

  return 0;
}

/* Select the elements in the upper triangle of the matrix. All other elements
 * will be set to zero.
 */
void
qdm_matrix_select_upper_triangle(gsl_matrix *m)
{
  for (size_t i = 0; i < m->size1; i++) {
    for (size_t j = 0; j < m->size2; j++) {
      if (j < i) {
        gsl_matrix_set(m, i, j, 0);
      }
    }
  }
}

int
qdm_matrix_to_csc_matrix(csc **result, gsl_matrix *input)
{
  int status = 0;

  gsl_spmatrix *sm = gsl_spmatrix_alloc(input->size1, input->size2);

  status = gsl_spmatrix_d2sp(sm, input);
  if (status != 0) {
    goto cleanup_sm;
  }

  gsl_spmatrix *sm_csc = gsl_spmatrix_compress(sm, GSL_SPMATRIX_CSC);

  int m = sm_csc->size1;
  int n = sm_csc->size2;
  int nzmax = gsl_spmatrix_nnz(sm_csc);

  size_t x_size = sizeof(double) * nzmax;
  double *x = malloc(x_size);
  memcpy(x, sm_csc->data, x_size);

  size_t i_size = sizeof(int) * nzmax;
  int *i = malloc(i_size);
  memcpy(i, sm_csc->i, i_size);

  size_t p_size = sizeof(int) * (n + 1);
  int *p = malloc(p_size);
  memcpy(p, sm_csc->p, p_size);

  *result = csc_matrix(
      m,     // m     First dimension (rows)
      n,     // n     Second dimension (columns)
      nzmax, // nzmax Maximum number of nonzero elements
      x,     // x     Vector of data (size nzmax)
      i,     // i     Vector of row indices (size nzmax)
      p      // p     Vector of column pointers (size n+1)
  );

  gsl_spmatrix_free(sm_csc);

cleanup_sm:
  gsl_spmatrix_free(sm);

  return status;
}

gsl_vector *
qdm_vector_copy(const gsl_vector *src)
{
  gsl_vector *dst = gsl_vector_alloc(src->size);
  gsl_vector_memcpy(dst, src);

  return dst;
}

gsl_matrix *
qdm_matrix_copy(const gsl_matrix *src)
{
  gsl_matrix *dst = gsl_matrix_alloc(src->size1, src->size2);
  gsl_matrix_memcpy(dst, src);

  return dst;
}

int create_hd5(
    const char *file_path,
    const char *group_path
)
{
  int status = 0;

  hid_t file  = -1;
  hid_t group = -1;
  hid_t gcpl  = -1;

  /* Save old error handler */
  herr_t (*old_func)(hid_t, void*) = NULL;
  void *old_client_data = NULL;
  H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  file = H5Fopen(file_path, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file < 0) {
    file = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  }
  if (file < 0) {
    status = file;
    goto cleanup;
  }

  group = H5Gopen(file, group_path, H5P_DEFAULT);
  if (group < 0) {
    gcpl = H5Pcreate(H5P_LINK_CREATE);
    if (gcpl < 0) {
      status = gcpl;
      goto cleanup;
    }

    status = H5Pset_create_intermediate_group(gcpl, 1);
    if (status < 0) {
      goto cleanup;
    }

    group = H5Gcreate(file, group_path, gcpl, H5P_DEFAULT, H5P_DEFAULT);

    if (group < 0) {
      status = group;
      goto cleanup;
    }
  }

cleanup:
  if (gcpl >= 0) {
    H5Pclose(gcpl);
  }

  if (group >= 0) {
    H5Gclose(group);
  }

  if (file >= 0) {
    H5Fclose(file);
  }

  /* Restore previous error handler */
  H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

  if (status < 0) {
    H5Eprint(H5E_DEFAULT, stderr);
  }

  return status;
}

int
qdm_vector_hd5_read(
    hid_t id,
    const char *name,
    gsl_vector **v
)
{
  int status = 0;

  int rank = 0;

  status = H5LTget_dataset_ndims(id, name, &rank);
  if (status < 0) {
    return status;
  }

  hsize_t dims[rank];

  status = H5LTget_dataset_info(id, name, dims, NULL, NULL);
  if (status < 0) {
    return status;
  }

  size_t size = 1;
  for (int i = 0; i < rank; i++) {
    size *= dims[i];
  }

  gsl_vector *tmp = gsl_vector_alloc(size);

  status = H5LTread_dataset_double(id, name, tmp->data);
  if (status < 0) {
    gsl_vector_free(tmp);

    return status;
  }

  *v = tmp;

  return status;
}

int
qdm_vector_hd5_write(
    hid_t id,
    const char *name,
    const gsl_vector *v
)
{
  int status = 0;

  hid_t datatype  = -1;
  hid_t dataspace = -1;
  hid_t dataset   = -1;
  hid_t dcpl      = -1;

  if (v == NULL) {
    goto cleanup;
  }

  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);

  status = H5Tset_order(datatype, H5T_ORDER_LE);
  if (status != 0) {
    goto cleanup;
  }

  hsize_t dims[1] = {
    v->size,
  };
  dataspace = H5Screate_simple(1, dims, NULL);
  if (dataspace < 0) {
    status = dataspace;
    goto cleanup;
  }

  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  if (dcpl < 0) {
    status = dcpl;

    goto cleanup;
  }

  /* Only enable compression on larger vectors. */
  /*
  if (v->size > 1024) {
    status = H5Pset_deflate(dcpl, 9);
    if (status != 0) {
      goto cleanup;
    }

    hsize_t chunk_dims[1] = {
      v->size,
    };
    status = H5Pset_chunk(dcpl, 1, chunk_dims);
    if (status != 0) {
      goto cleanup;
    }
  }
  */

  dataset = H5Dcreate(id, name, datatype, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  if (dataset < 0) {
    status = dataset;
    goto cleanup;
  }

  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v->data);
  if (status != 0) {
    goto cleanup;
  }

cleanup: 
  if (dcpl >= 0) {
    H5Pclose(dcpl);
  }

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
qdm_matrix_hd5_read(
    hid_t id,
    const char *name,
    gsl_matrix **m
)
{
  int status = 0;

  int rank = 0;

  status = H5LTget_dataset_ndims(id, name, &rank);
  if (status < 0) {
    return status;
  }

  if (rank < 2) {
    return -1;
  }

  hsize_t dims[rank];

  status = H5LTget_dataset_info(id, name, dims, NULL, NULL);
  if (status < 0) {
    return status;
  }

  size_t size1 = dims[0];
  size_t size2 = 1;
  for (int i = 1; i < rank; i++) {
    size2 *= dims[i];
  }

  gsl_matrix *tmp = gsl_matrix_alloc(size1, size2);

  status = H5LTread_dataset_double(id, name, tmp->data);
  if (status < 0) {
    gsl_matrix_free(tmp);

    return status;
  }

  *m = tmp;

  return status;
}

int
qdm_matrix_hd5_write(
    hid_t id,
    const char *name,
    const gsl_matrix *m
)
{
  int status = 0;

  hid_t datatype  = -1;
  hid_t dataspace = -1;
  hid_t dataset   = -1;
  hid_t dcpl      = -1;

  if (m == NULL) {
    goto cleanup;
  }

  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);

  status = H5Tset_order(datatype, H5T_ORDER_LE);
  if (status != 0) {
    goto cleanup;
  }

  hsize_t dims[2] = {
    m->size1,
    m->size2,
  };
  dataspace = H5Screate_simple(2, dims, NULL);
  if (dataspace < 0) {
    status = dataspace;

    goto cleanup;
  }

  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  if (dcpl < 0) {
    status = dcpl;

    goto cleanup;
  }

  /* Only enable compression on larger matrices. */
  /*
  if (m->size1 * m->size2 > 1024) {
    status = H5Pset_deflate(dcpl, 9);
    if (status != 0) {
      goto cleanup;
    }

    hsize_t chunk_dims[2] = {
      m->size1,
      m->size2
    };
    status = H5Pset_chunk(dcpl, 2, chunk_dims);
    if (status != 0) {
      goto cleanup;
    }
  }
  */

  dataset = H5Dcreate(id, name, datatype, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  if (dataset < 0) {
    status = dataset;

    goto cleanup;
  }

  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, m->data);
  if (status != 0) {
    goto cleanup;
  }

cleanup:
  if (dcpl >= 0) {
    H5Pclose(dcpl);
  }

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

gsl_vector *
qdm_matrix_filter(
    const gsl_matrix *m,
    size_t needle_column,
    double needle_value,
    size_t select_column
)
{
  size_t size = 0;

  for (size_t i = 0; i < m->size1; i++) {
    if (gsl_matrix_get(m, i, needle_column) == needle_value) {
      size++;
    }
  }

  gsl_vector *v = gsl_vector_alloc(size);

  size_t j = 0;
  for (size_t i = 0; i < m->size1; i++) {
    double select_value = gsl_matrix_get(m, i, select_column);

    if (gsl_matrix_get(m, i, needle_column) == needle_value) {
      gsl_vector_set(v, j, select_value);
      j++;
    }
  }

  return v; 
}
