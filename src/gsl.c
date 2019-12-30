/* Utility functions for working with GSL data types. */

#include <math.h>
#include <stdio.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_statistics.h>
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
qdm_vector_csv_fwrite(FILE *f, gsl_vector *v)
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
qdm_matrix_csv_fwrite(FILE *f, gsl_matrix *m)
{
  for (size_t i = 0; i < m->size1; i++) {
    gsl_vector_view row = gsl_matrix_row(m, i);
    qdm_vector_csv_fwrite(f, &row.vector);
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

gsl_vector *
qdm_vector_quantile(gsl_vector *data, gsl_vector *probs)
{
  gsl_vector *quantiles = gsl_vector_alloc(probs->size);

  for (size_t i = 0; i < probs->size; i++) {
    double q = gsl_stats_quantile_from_sorted_data(
        data->data,
        data->stride,
        data->size,
        gsl_vector_get(probs, i)
    );

    gsl_vector_set(quantiles, i, q);
  }

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
