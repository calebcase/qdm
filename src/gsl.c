/* Utility functions for working with GSL data types. */

#include <math.h>
#include <stdio.h>

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

void
qdm_matrix_csv_fwrite(FILE *f, gsl_matrix *m)
{
  for (size_t i = 0; i < m->size1; i++) {
    gsl_vector_view row = gsl_matrix_row(m, i);
    qdm_vector_csv_fwrite(f, &row.vector);
  }
}
