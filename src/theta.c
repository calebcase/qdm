#include <float.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
#include <osqp/osqp.h>

#include "qdm.h"

#define QDM_THETA_MIN_THRESHOLD 0.0001

int
qdm_theta_optimize(
    gsl_vector *result,
    gsl_vector *emperical_quantiles,
    gsl_matrix *ix
)
{
  int status = 0;

  gsl_matrix *p = gsl_matrix_alloc(ix->size2, ix->size2);
  gsl_matrix *a = gsl_matrix_alloc(ix->size2, ix->size2);
  gsl_vector *q = gsl_vector_alloc(ix->size2);
  gsl_vector *l = gsl_vector_alloc(ix->size2);
  gsl_vector *u = gsl_vector_alloc(ix->size2);

  // FIXME check c_malloc...
  OSQPData *data         = c_malloc(sizeof(OSQPData));
  OSQPSettings *settings = c_malloc(sizeof(OSQPSettings));
  OSQPWorkspace *work    = NULL;

  status = qdm_matrix_tmm(ix, p);
  if (status != 0) {
    goto cleanup;
  }

  qdm_matrix_select_upper_triangle(p);

  gsl_matrix_set_identity(a);

  gsl_vector_set_all(l, QDM_THETA_MIN_THRESHOLD);
  gsl_vector_set_all(u, INFINITY);

  status = gsl_blas_dgemv(CblasTrans, -1.0, ix, emperical_quantiles, 0, q);
  if (status != 0) {
    goto cleanup;
  }

  /* Solve... */
  data->n = ix->size2;
  data->m = ix->size2;

  status = qdm_matrix_to_csc_matrix(&data->P, p);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_matrix_to_csc_matrix(&data->A, a);
  if (status != 0) {
    goto cleanup;
  }

  data->q = q->data;
  data->l = l->data;
  data->u = u->data;

  osqp_set_default_settings(settings);
  settings->verbose = 0;

  status = osqp_setup(&work, data, settings);
  if (status != 0) {
    goto cleanup;
  }

  status = osqp_solve(work);
  if (status != 0) {
    goto cleanup;
  }

  gsl_vector_view beta0 = gsl_vector_view_array(work->solution->x, ix->size2);
  status = gsl_vector_memcpy(result, &beta0.vector);
  if (status != 0) {
    goto cleanup;
  }

cleanup:
  if (data != NULL) {
    c_free(data->A);
    c_free(data->P);
    c_free(data);
  }

  c_free(settings);
  osqp_cleanup(work);

  gsl_vector_free(u);
  gsl_vector_free(l);
  gsl_vector_free(q);
  gsl_matrix_free(a);
  gsl_matrix_free(p);

  return status;
}

void
qdm_theta_matrix_constrain(
    gsl_matrix *theta,
    double min
)
{
  /* First row is constrained to be greater than zero. */
  for (size_t j = 0; j < theta->size2; j++) {
    if (gsl_matrix_get(theta, 0, j) < 0) {
      gsl_matrix_set(theta, 0, j, min);
    }
  }

  /* Remaining cells are set to zero if the sum of the column is negative.
   *
   * NOTE: The first column is not modified.
   */
  for (size_t i = 1; i < theta->size1; i++) {
    for (size_t j = 1; j < theta->size2; j++) {
      gsl_vector_view column = gsl_matrix_column(theta, j);
      double sum = qdm_vector_sum(&column.vector);

      if (sum < 0) {
        gsl_matrix_set(theta, i, j, 0);
      }
    }
  }
}
