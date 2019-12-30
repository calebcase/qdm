#include <math.h>
#include <stdio.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_vector.h>

#include <argtable2.h>

#include "qdm.h"

void
qdm_ispline_vector(
    gsl_vector *result,
    const double tau,
    const size_t spline_df,
    const gsl_vector *knots
)
{
  size_t bin = qdm_vector_search(knots, tau);

  for (size_t m = 0; m < result->size - 1; m++) {
    double v;

    if (bin < m) {
      v = 0;
    } else if (bin - spline_df + 1 > m) {
      v = 1;
    } else if (bin == m) {
      double n = gsl_sf_pow_int(tau - gsl_vector_get(knots, m), 3);
      double d1 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m);
      double d2 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m);
      double d3 = gsl_vector_get(knots, m + 1) - gsl_vector_get(knots, m);

      v = n / (d1 * d2 * d3);
    } else if (bin == m + 1) {
      double i1n = 3 * tau;
      double i1d = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m + 1);
      double i1 = i1n / i1d;

      double i2n = -gsl_sf_pow_int(tau - gsl_vector_get(knots, m), 3);
      double i2d1 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m);
      double i2d2 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m);
      double i2d3 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m + 1);
      double i2 = i2n / (i2d1 * i2d2 * i2d3);

      double i3n = gsl_sf_pow_int(gsl_vector_get(knots, m + 3) - tau, 3);
      double i3d1 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m + 1);
      double i3d2 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m);
      double i3d3 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m + 1);
      double i3 = i3n / (i3d1 * i3d2 * i3d3);

      double i4an = -3 * gsl_vector_get(knots, m + 1);
      double i4ad = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m + 1);
      double i4a = i4an / i4ad;
      double i4bn = gsl_sf_pow_int(gsl_vector_get(knots, m + 1) - gsl_vector_get(knots, m), 3);
      double i4bd1 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m);
      double i4bd2 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m);
      double i4bd3 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m + 1);
      double i4b = i4bn / (i4bd1 * i4bd2 * i4bd3);
      double i4cn = gsl_sf_pow_int(gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m + 1), 3);
      double i4cd1 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m + 1);
      double i4cd2 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m);
      double i4cd3 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m + 1);
      double i4c = i4cn / (i4cd1 * i4cd2 * i4cd3);
      double i4 = i4a + i4b - i4c;

      double i5 = 0;
      if (m != 1) {
        double i5n = gsl_sf_pow_int(gsl_vector_get(knots, m + 1) - gsl_vector_get(knots, m), 3);
        double i5d1 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m);
        double i5d2 = gsl_vector_get(knots, m + 2) - gsl_vector_get(knots, m);
        double i5d3 = gsl_vector_get(knots, m + 1) - gsl_vector_get(knots, m);
        i5 = i5n / (i5d1 * i5d2 * i5d3);
      }

      v = i1 + i2 + i3 + i4 + i5;
    } else {
      double n = gsl_sf_pow_int(gsl_vector_get(knots, m + 3) - tau, 3);
      double d1 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m);
      double d2 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m + 1);
      double d3 = gsl_vector_get(knots, m + 3) - gsl_vector_get(knots, m + 2);
      double d = d1 * d2 * d3;

      v = 1 - n / d;
    }

    gsl_vector_set(result, 0, 1);
    gsl_vector_set(result, m + 1, v);
  }
}

void
qdm_ispline_matrix(
    gsl_matrix *result,
    const gsl_vector *taus,
    const size_t spline_df,
    const gsl_vector *knots
)
{
  for (size_t i = 0; i < result->size1; i++) {
    gsl_vector_view row = gsl_matrix_row(result, i);
    qdm_ispline_vector(&row.vector, gsl_vector_get(taus, i), spline_df, knots);
  }
}
