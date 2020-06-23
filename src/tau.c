#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "qdm.h"

#define TAU_RESET_L_SIZE 11
#define TAU_RESET_M_SIZE  8
#define TAU_RESET_H_SIZE 11
#define TAU_RESET_SIZE (TAU_RESET_L_SIZE + TAU_RESET_M_SIZE + TAU_RESET_H_SIZE)

#define TAU_RESET_STEP 0.0001
#define TAU_RESET_GAP  0.001

#define TAU_TABLE_SIZE 1048576

#define TAU_INITIAL_GUESS 0.5
#define TAU_EPS 0.001
#define TAU_ITERATION_MAX 1000

static
void
qdm_tau_reset_setup(
    qdm_tau *t
)
{
  t->reset = gsl_vector_alloc(TAU_RESET_SIZE);

  size_t offset = 0;

  gsl_vector_view view;

  view = gsl_vector_subvector(t->reset, offset, TAU_RESET_L_SIZE);
  qdm_vector_set_seq(&view.vector, t->low, t->low + (TAU_RESET_L_SIZE - 1) * TAU_RESET_STEP);
  offset += TAU_RESET_L_SIZE;

  view = gsl_vector_subvector(t->reset, offset, TAU_RESET_M_SIZE);
  qdm_vector_set_seq(&view.vector, t->low + TAU_RESET_GAP, t->high - TAU_RESET_GAP);
  offset += TAU_RESET_M_SIZE;

  view = gsl_vector_subvector(t->reset, offset, TAU_RESET_H_SIZE);
  qdm_vector_set_seq(&view.vector, t->high - (TAU_RESET_H_SIZE - 1) * TAU_RESET_STEP, t->high);
  offset += TAU_RESET_H_SIZE;
}

static
void
qdm_tau_table_setup(
    qdm_tau *t
)
{
  if (!t->use_table) {
    t->ispline_table = NULL;
    t->mspline_table = NULL;

    return;
  }

  t->ispline_table = gsl_matrix_alloc(TAU_TABLE_SIZE, (t->knots->size - t->spline_df) + 1);
  t->mspline_table = gsl_matrix_alloc(TAU_TABLE_SIZE, (t->knots->size - t->spline_df) + 1);

  double tau = 0;
  gsl_vector_view row;
  for (size_t i = 0; i < TAU_TABLE_SIZE; i++) {
    tau = (double)(i) / TAU_TABLE_SIZE;

    row = gsl_matrix_row(t->ispline_table, i);
    qdm_ispline_vector(
        &row.vector,
        tau,
        t->spline_df,
        t->knots
    );

    row = gsl_matrix_row(t->mspline_table, i);
    qdm_mspline_vector(
        &row.vector,
        tau,
        t->spline_df,
        t->knots
    );
  }
}

qdm_tau *
qdm_tau_alloc(
    int use_table,

    double low,
    double high,

    size_t spline_df,
    const gsl_vector *knots
)
{
  qdm_tau *t = malloc(sizeof(qdm_tau));

  if (use_table == 1) {
    t->use_table = true;
  } else {
    t->use_table = false;
  }

  t->low = low;
  t->high = high;

  t->spline_df = spline_df;
  t->knots = qdm_vector_copy(knots);

  qdm_tau_reset_setup(t);
  qdm_tau_table_setup(t);

  return t;
}

void
qdm_tau_free(
    qdm_tau *t
)
{
  if (t == NULL) {
    return;
  }

  t->use_table = 0;

  t->low = 0;
  t->high = 0;
  
  t->spline_df = 0;
  t->knots = NULL;

  gsl_vector_free(t->reset);
  t->reset = NULL;

  gsl_matrix_free(t->ispline_table);
  t->ispline_table = NULL;

  gsl_matrix_free(t->mspline_table);
  t->mspline_table = NULL;

  free(t);
}

double
qdm_tau_ispline_mmm(
    const qdm_tau *t,
    const double value,
    const gsl_vector *mmm
)
{
  double result = 0;

  if (t->use_table) {
    size_t i = (size_t)fmin(floor(value * TAU_TABLE_SIZE), TAU_TABLE_SIZE - 1);
    gsl_vector_view ispline = gsl_matrix_row(t->ispline_table, i);

    gsl_blas_ddot(&ispline.vector, mmm, &result);
  } else {
    size_t m = (t->knots->size - t->spline_df) + 1;
    double ispline_data[m];
    gsl_vector_view ispline = gsl_vector_view_array(ispline_data, m);

    qdm_ispline_vector(&ispline.vector, value, t->spline_df, t->knots);
    gsl_blas_ddot(&ispline.vector, mmm, &result);
  }

  return result;
}

double
qdm_tau_mspline_mmm(
    const qdm_tau *t,
    const double value,
    const gsl_vector *mmm
)
{
  double result = 0;

  if (t->use_table) {
    size_t i = (size_t)fmin(floor(value * TAU_TABLE_SIZE), TAU_TABLE_SIZE - 1);
    gsl_vector_view mspline = gsl_matrix_row(t->mspline_table, i);

    gsl_blas_ddot(&mspline.vector, mmm, &result);
  } else {
    size_t m = (t->knots->size - t->spline_df) + 1;
    double mspline_data[m];
    gsl_vector_view mspline = gsl_vector_view_array(mspline_data, m);

    qdm_mspline_vector(&mspline.vector, value, t->spline_df, t->knots);
    gsl_blas_ddot(&mspline.vector, mmm, &result);
  }

  return result;
}

double
qdm_tau_find(
    const qdm_tau *t,
    const double v,
    const gsl_vector *mmm
)
{
  double tau = TAU_INITIAL_GUESS;
  double qi_u = 0;
  double qm_u = 0;

  size_t reset_i = 0;
  for (size_t i = 0; i < TAU_ITERATION_MAX; i++) {
    if (reset_i >= t->reset->size) {
      tau = TAU_INITIAL_GUESS;

      break;
    }

    /* Calculate position... */
    qi_u = qdm_tau_ispline_mmm(t, tau, mmm);

    /* Calculate slope... */
    qm_u = qdm_tau_mspline_mmm(t, tau, mmm);

    /* Check if the update is within our desired interval [0, 1]. */
    double update = tau - (qi_u - v) / qm_u;
    if (update < 0 || update > 1) {
      tau = gsl_vector_get(t->reset, reset_i);
      reset_i++;

      continue;
    } else {
      tau = update;
    }

    if (fabs(qi_u - v) <= TAU_EPS) {
      break;
    }
  }

  return tau;
}
