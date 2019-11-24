#include <math.h>
#include <stdio.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_vector.h>

gsl_vector *
gsl_vector_seq(double from, double to, double by)
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

gsl_vector *
gsl_vector_M_knots(size_t spline_df, gsl_vector *knots_inter)
{
  size_t size = spline_df * 2 + knots_inter->size;
  gsl_vector *M_knots = gsl_vector_alloc(size);

  // Fill front with zeros.
  for (size_t i = 0; i < spline_df; i++) {
    gsl_vector_set(M_knots, i, 0);
  }

  // Fill middle with knots.
  gsl_vector_view view = gsl_vector_subvector(M_knots, spline_df, knots_inter->size);
  gsl_blas_dcopy(knots_inter, &view.vector);

  // Fill end with ones.
  for (size_t i = spline_df + knots_inter->size; i < M_knots->size; i++) {
    gsl_vector_set(M_knots, i, 1);
  }

  return M_knots;
}

size_t
gsl_vector_search(const gsl_vector *v, double needle)
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
gsl_vector_mspline(gsl_vector *result, const double tau, const size_t spline_df, const gsl_vector *knots)
{
  size_t bin = gsl_vector_search(knots, tau);

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
gsl_matrix_mspline(gsl_matrix *result, const size_t spline_df, const gsl_vector *x, const gsl_vector *knots)
{
  for (size_t i = 0; i < result->size1; i++) {
    gsl_vector_view row = gsl_matrix_row(result, i);
    gsl_vector_mspline(&row.vector, gsl_vector_get(x, i), spline_df, knots);
  }
}

void
gsl_vector_ispline(gsl_vector *result, const double tau, const size_t spline_df, const gsl_vector *knots)
{
  size_t bin = gsl_vector_search(knots, tau);

  for (size_t m = 0; m < result->size; m++) {
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

    gsl_vector_set(result, m, v);
  }
}

void
gsl_matrix_ispline(gsl_matrix *result, const size_t spline_df, const gsl_vector *x, const gsl_vector *knots)
{
  for (size_t i = 0; i < result->size1; i++) {
    gsl_vector_view row = gsl_matrix_row(result, i);
    gsl_vector_ispline(&row.vector, gsl_vector_get(x, i), spline_df, knots);
  }
}

void
gsl_vector_csv_fwrite(FILE *f, gsl_vector *v)
{
  for (size_t i = 0; i < v->size; i++) {
    fprintf(f, "%.17g", gsl_vector_get(v, i));
    if (i < v->size - 1) {
      fprintf(f, ",");
    }
  }
  fprintf(f, "\n");
}

void
gsl_matrix_csv_fwrite(FILE *f, gsl_matrix *m)
{
  for (size_t i = 0; i < m->size1; i++) {
    gsl_vector_view row = gsl_matrix_row(m, i);
    gsl_vector_csv_fwrite(f, &row.vector);
  }
}

int
main()
{
  size_t spline_df = 3;

  printf("knots_inter\n");
  gsl_vector *knots_inter = gsl_vector_seq(0.1, 0.9, 0.2);
  gsl_vector_csv_fwrite(stdout, knots_inter);

  printf("x\n");
  gsl_vector *x = gsl_vector_seq(0.01, 0.99, 0.01);
  gsl_vector_csv_fwrite(stdout, x);

  printf("knots\n");
  gsl_vector *knots = gsl_vector_M_knots(spline_df, knots_inter);
  gsl_vector_csv_fwrite(stdout, knots);

  size_t M = spline_df + knots_inter->size;
  printf("M\n%zi\n", M);

  gsl_matrix *MX = gsl_matrix_alloc(x->size, M);
  gsl_matrix_mspline(MX, spline_df, x, knots);

  printf("MX\n");
  gsl_matrix_csv_fwrite(stdout, MX);

  gsl_matrix *IX = gsl_matrix_alloc(x->size, M);
  gsl_matrix_ispline(IX, spline_df, x, knots);

  printf("IX\n");
  gsl_matrix_csv_fwrite(stdout, IX);

  return 0;
}
