#include <math.h>
#include <stdio.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <qdm.h>

gsl_vector *
qdm_knots_vector(size_t spline_df, gsl_vector *interior_knots)
{
  size_t size = spline_df * 2 + interior_knots->size;
  gsl_vector *knots = gsl_vector_alloc(size);

  // Fill front with zeros.
  for (size_t i = 0; i < spline_df; i++) {
    gsl_vector_set(knots, i, 0);
  }

  // Fill middle with interior knots.
  gsl_vector_view view = gsl_vector_subvector(knots, spline_df, interior_knots->size);
  gsl_blas_dcopy(interior_knots, &view.vector);

  // Fill end with ones.
  for (size_t i = spline_df + interior_knots->size; i < knots->size; i++) {
    gsl_vector_set(knots, i, 1);
  }

  return knots;
}

int
qdm_knots_rss(
  double *rss,

  gsl_vector *sorted_data,
  gsl_vector *middle,
  gsl_vector *interior_knots,

  size_t spline_df
)
{
  int status = 0;

  gsl_vector *m_knots = qdm_knots_vector(spline_df, interior_knots);

  size_t m = spline_df + interior_knots->size;
  gsl_matrix *ix = gsl_matrix_alloc(middle->size, m + 1);

  gsl_vector *emperical_quantiles = qdm_vector_quantile(sorted_data, middle);
  gsl_vector *theta = gsl_vector_alloc(ix->size2);
  gsl_vector *spline_quantiles = gsl_vector_alloc(emperical_quantiles->size);

  /* Build the true quantile function. */
  qdm_ispline_matrix(ix, middle, spline_df, m_knots);

  /* Check if the matrix is full rank. */
  double det = 0;
  status = qdm_matrix_det_tmm(ix, &det);
  if (status != 0) {
    goto cleanup;
  }

  if (det <= 0) {
    goto cleanup;
  }

  status = qdm_theta_optimize(theta, emperical_quantiles, ix);
  if (status != 0) {
    goto cleanup;
  }

  status = gsl_blas_dgemv(CblasNoTrans , 1.0, ix, theta, 0.0, spline_quantiles);
  if (status != 0) {
    goto cleanup;
  }

  *rss = qdm_vector_rss(emperical_quantiles, spline_quantiles);

cleanup:
  gsl_vector_free(spline_quantiles);
  gsl_vector_free(theta);
  gsl_vector_free(emperical_quantiles);
  gsl_matrix_free(ix);
  gsl_vector_free(m_knots);

  return status;
}

int
qdm_knots_optimize(
  gsl_vector *result,

  gsl_vector *sorted_data,
  gsl_vector *middle,
  gsl_vector *possible_knots,

  size_t iterate_n,
  size_t spline_df
)
{
  int status = 0;
  double rss = INFINITY;
  double rss_found = 0;

  gsl_vector *interior_knots = gsl_vector_alloc(result->size);

  for (size_t i = 0; i < iterate_n; i++) {
    status = gsl_ran_choose(
        qdm_state->rng,
        interior_knots->data,
        interior_knots->size,
        possible_knots->data,
        possible_knots->size,
        sizeof(double)
    );
    if (status != 0) {
      goto cleanup;
    }

    status = qdm_knots_rss(
        &rss_found,

        sorted_data,
        middle,
        interior_knots,

        spline_df
    );
    if (status != 0) {
      goto cleanup;
    }

    if (rss_found < rss) {
      rss = rss_found;

      status = gsl_vector_memcpy(result, interior_knots);
      if (status != 0) {
        goto cleanup;
      }
    }
  }

cleanup:
  gsl_vector_free(interior_knots);

  return status;
}
