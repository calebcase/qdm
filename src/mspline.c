#include <math.h>
#include <stdio.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_vector.h>

#include <argtable2.h>

#include "qdm.h"

void
qdm_mspline_vector(gsl_vector *result, const double tau, const size_t spline_df, const gsl_vector *knots)
{
  // FIXME First column needs to be 0
  size_t bin = qdm_vector_search(knots, tau);

  for (size_t m = 0; m < result->size; m++) {
    double v;

    if (bin < m) {
      v = 0;
    } else if (bin - spline_df + 1 > m) {
      v = 0;
    } else if (bin == m) {
      double n = 3 * gsl_sf_pow_int(tau - gsl_vector_get(knots, m), 2);
      double d1 = gsl_vector_get(knots, m + 1) - gsl_vector_get(knots, m);
      double d2 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m);
      double d3 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m);

      v = n / (d1 * d2 * d3);
    } else if (bin == m + 1) {
      double i1 = 1 / (gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m + 1));

      double i2n = gsl_sf_pow_int(tau - gsl_vector_get(knots, m), 2);
      double i2d1 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m);
      double i2d2 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m);
      double i2d3 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m + 1);
      double i2 = i2n / (i2d1 * i2d2 * i2d3);

      double i3n = gsl_sf_pow_int(gsl_vector_get(knots, m + 3) - tau, 2);
      double i3d1 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m);
      double i3d2 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m + 1);
      double i3d3 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m + 1);
      double i3 = i3n / (i3d1 * i3d2 * i3d3);

      v = 3 * (i1 - i2 - i3);
    } else {
      double n = 3 * gsl_sf_pow_int(gsl_vector_get(knots, m + 3) - tau, 2);
      double d1 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m);
      double d2 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m + 1);
      double d3 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m + 2);
      double d = d1 * d2 * d3;

      v = n / d;
    }

    gsl_vector_set(result, m, v);
  }
}

void
qdm_mspline_matrix(gsl_matrix *result, const size_t spline_df, const gsl_vector *x, const gsl_vector *knots)
{
  for (size_t i = 0; i < result->size1; i++) {
    gsl_vector_view row = gsl_matrix_row(result, i);
    qdm_mspline_vector(&row.vector, gsl_vector_get(x, i), spline_df, knots);
  }
}
