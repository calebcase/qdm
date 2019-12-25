#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

gsl_vector *
qdm_knots_vector(size_t spline_df, gsl_vector *knots_inter)
{
  size_t size = spline_df * 2 + knots_inter->size;
  gsl_vector *knots = gsl_vector_alloc(size);

  // Fill front with zeros.
  for (size_t i = 0; i < spline_df; i++) {
    gsl_vector_set(knots, i, 0);
  }

  // Fill middle with knots.
  gsl_vector_view view = gsl_vector_subvector(knots, spline_df, knots_inter->size);
  gsl_blas_dcopy(knots_inter, &view.vector);

  // Fill end with ones.
  for (size_t i = spline_df + knots_inter->size; i < knots->size; i++) {
    gsl_vector_set(knots, i, 1);
  }

  return knots;
}
