#define _GNU_SOURCE
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <argtable2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_minmax.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics_double.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <qdm.h>

#define MAX_ARG_ERRORS 20

#define FUTURE_YEAR  0
#define FUTURE_MONTH 1
#define FUTURE_DAY   2
#define FUTURE_VALUE 3

#define HISTORICAL_YEAR     0
#define HISTORICAL_MONTH    1
#define HISTORICAL_DAY      2
#define HISTORICAL_OBSERVED 3
#define HISTORICAL_VALUE    4

int
qdm_double_write(
    hid_t id,
    const char *name,
    double value
)
{
  hsize_t dims[1] = {1};
  double data[1] = {value};

  return H5LTmake_dataset_double(id, name, 1, dims, data);
}

int
qdm_int_write(
    hid_t id,
    const char *name,
    int value
)
{
  hsize_t dims[1] = {1};
  int data[1] = {value};

  return H5LTmake_dataset_int(id, name, 1, dims, data);
}

typedef struct {
  int rng_seed;

  int acc_check;
  int burn;
  int iter;
  int knot_try;
  int spline_df;
  int thin;

  double month;
  double knot;

  double tau_high;
  double tau_low;

  double theta_min;
  double theta_tune_sd;

  double xi_high;
  double xi_low;
  double xi_prior_mean;
  double xi_prior_var;
  double xi_tune_sd;
} qdm_parameters;

void
qdm_parameters_fprint(
    FILE *f,
    const qdm_parameters *p
)
{
  const char *prefix = "    ";

  fprintf(f, "%smonth: %f\n", prefix, p->month);
  fprintf(f, "%sknot: %f\n",  prefix, p->knot);

  fprintf(f, "%stau_high: %f\n", prefix, p->tau_high);
  fprintf(f, "%stau_low: %f\n",  prefix, p->tau_low);
}

int
qdm_parameters_write(
    hid_t id,
    const qdm_parameters *p
)
{
  int status = 0;

#define WRITE_INT(n) \
  status = qdm_int_write(id, #n, p->n); \
  if (status != 0) { \
    return status; \
  }

#define WRITE_DOUBLE(n) \
  status = qdm_double_write(id, #n, p->n); \
  if (status != 0) { \
    return status; \
  }

  WRITE_INT(rng_seed);

  WRITE_INT(acc_check);
  WRITE_INT(burn);
  WRITE_INT(iter);
  WRITE_INT(knot_try);
  WRITE_INT(spline_df);
  WRITE_INT(thin);

  WRITE_DOUBLE(month);
  WRITE_DOUBLE(knot);

  WRITE_DOUBLE(tau_high);
  WRITE_DOUBLE(tau_low);

  WRITE_DOUBLE(theta_min);
  WRITE_DOUBLE(theta_tune_sd);

  WRITE_DOUBLE(xi_high);
  WRITE_DOUBLE(xi_low);
  WRITE_DOUBLE(xi_prior_mean);
  WRITE_DOUBLE(xi_prior_var);
  WRITE_DOUBLE(xi_tune_sd);

#undef WRITE_DOUBLE
#undef WRITE_INT

  return status;
}

typedef struct {
  // Parameters
  qdm_parameters parameters;

  // Years (pre-scaling)
  double years_min;
  double years_max;

  // Intermediate Results
  size_t intermediates_count;
  qdm_intermediate_result *intermediates;

  // Diagnostics
  double waic;
  double pwaic;
  double dic;
  double pd;

  // Results
  gsl_matrix *theta_bar;
  double xi_low_bar;
  double xi_high_bar;

  gsl_matrix *theta_acc;
  gsl_vector *xi_acc;

  gsl_vector *m_knots;

  gsl_matrix *bias_corrected;

  // Performance Metrics
  double elapsed;
} qdm_evaluation;

qdm_evaluation *
qdm_evaluation_new(const qdm_parameters *parameters)
{
  qdm_evaluation *e = malloc(sizeof(qdm_evaluation));

  e->parameters = *parameters;

  e->years_min = 0;
  e->years_max = 0;

  e->intermediates_count = parameters->iter / parameters->thin - 1;
  e->intermediates = calloc(sizeof(qdm_intermediate_result), e->intermediates_count);
 
  e->waic = 0;
  e->pwaic = 0;
  e->dic = 0;
  e->pd = 0;

  e->theta_bar = NULL;
  e->xi_low_bar = 0;
  e->xi_high_bar = 0;

  e->theta_acc = NULL;
  e->xi_acc = NULL;

  e->m_knots = NULL;

  e->bias_corrected = NULL;

  e->elapsed = 0;

  return e;
}

void
qdm_evaluation_free(qdm_evaluation *e)
{
  gsl_matrix_free(e->bias_corrected);
  e->bias_corrected = NULL;

  gsl_vector_free(e->m_knots);
  e->m_knots = NULL;

  gsl_vector_free(e->xi_acc);
  e->xi_acc = NULL;

  gsl_matrix_free(e->theta_acc);
  e->theta_acc = NULL;

  gsl_matrix_free(e->theta_bar);
  e->theta_bar = NULL;

  for (size_t i = 0; i < e->intermediates_count; i++) {
    gsl_matrix_free(e->intermediates[i].theta);
    gsl_matrix_free(e->intermediates[i].theta_star);

    gsl_vector_free(e->intermediates[i].ll);
    gsl_vector_free(e->intermediates[i].tau);

    gsl_vector_free(e->intermediates[i].xi);
  }

  free(e->intermediates);
  e->intermediates = NULL;

  e->intermediates_count = 0;

  free(e);
}

void
qdm_evaluation_fprint(
    FILE *f,
    const qdm_evaluation *e
)
{
  const char *prefix = "  ";

  fprintf(f, "%sparameters:\n", prefix);
  qdm_parameters_fprint(f, &e->parameters);

  fprintf(f, "%syears_min: %f\n", prefix, e->years_min);
  fprintf(f, "%syears_max: %f\n", prefix, e->years_max);

  fprintf(f, "%sintermediates_count: %zu\n", prefix, e->intermediates_count);

  fprintf(f, "%swaic: %f\n", prefix, e->waic);
  fprintf(f, "%spwaic: %f\n", prefix, e->pwaic);
  fprintf(f, "%sdic: %f\n", prefix, e->dic);
  fprintf(f, "%spd: %f\n", prefix, e->pd);

  fprintf(f, "%stheta_bar:\n", prefix);
  qdm_matrix_csv_fwrite(f, e->theta_bar);

  fprintf(f, "%sxi_low_bar: %f\n", prefix, e->xi_low_bar);
  fprintf(f, "%sxi_high_bar: %f\n", prefix, e->xi_high_bar);

  fprintf(f, "%stheta_acc:\n", prefix);
  qdm_matrix_csv_fwrite(f, e->theta_acc);

  fprintf(f, "%sxi_acc:\n", prefix);
  qdm_vector_csv_fwrite(f, e->xi_acc);

  fprintf(f, "%sm_knots:\n", prefix);
  qdm_vector_csv_fwrite(f, e->m_knots);

  fprintf(f, "%selapsed: %f\n", prefix, e->elapsed);
}

int
qdm_evaluation_write(
    hid_t id,
    const qdm_evaluation *e
)
{
  int status = 0;

  hid_t parameters_group = -1;

#define WRITE_DOUBLE(n) \
  status = qdm_double_write(id, #n, e->n); \
  if (status != 0) { \
    goto cleanup; \
  }

  parameters_group = qdm_data_create_group(id, "parameters");
  if (parameters_group < 0) {
    status = parameters_group;

    goto cleanup;
  }

  status = qdm_parameters_write(parameters_group, &e->parameters);
  if (status != 0) {
    goto cleanup;
  }

  WRITE_DOUBLE(years_min);
  WRITE_DOUBLE(years_max);

  /*
  for (size_t i = 0; i < e->intermediates_count; i++) {
    status = qdm_data_intermediate_result_write(id, i, &e->intermediates[i]);
    if (status != 0) {
      goto cleanup;
    }
  }
  */
 
  WRITE_DOUBLE(waic);
  WRITE_DOUBLE(pwaic);
  WRITE_DOUBLE(dic);
  WRITE_DOUBLE(pd);

  status = qdm_matrix_hd5_write(id, "theta_bar", e->theta_bar);
  if (status != 0) {
    goto cleanup;
  }

  WRITE_DOUBLE(xi_high_bar);
  WRITE_DOUBLE(xi_low_bar);

  status = qdm_matrix_hd5_write(id, "theta_acc", e->theta_acc);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_vector_hd5_write(id, "xi_acc", e->xi_acc);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_vector_hd5_write(id, "m_knots", e->m_knots);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_matrix_hd5_write(id, "bias_corrected", e->bias_corrected);
  if (status != 0) {
    goto cleanup;
  }

  WRITE_DOUBLE(elapsed);

#undef WRITE_DOUBLE

cleanup:
  if (parameters_group >= 0) {
    H5Gclose(parameters_group);
  }

  return status;
}

int
qdm_evaluation_run(
    qdm_evaluation *e,
    const gsl_matrix *historical,
    const gsl_matrix *future
)
{
  int status = 0;

  clock_t started_at = clock();

  gsl_vector *years = NULL;
  gsl_vector *years_norm = NULL;
  gsl_vector *days = NULL;

  gsl_vector *values = NULL;
  gsl_vector *values_sorted = NULL;

  gsl_vector *interior_knots = NULL;

  gsl_vector *m_knots = NULL;
  gsl_matrix *theta = NULL;
  gsl_vector *ll = NULL;
  gsl_vector *tau = NULL;

  gsl_matrix *theta_tune = NULL;
  gsl_matrix *theta_acc = NULL;
  gsl_vector *xi_tune = NULL;
  gsl_vector *xi_acc = NULL;

  /* Initialize RNG */

  gsl_rng_set(qdm_state->rng, e->parameters.rng_seed);

  /* Normalize Years */

  {
    years = qdm_matrix_filter(future, FUTURE_MONTH, e->parameters.month, FUTURE_YEAR);
    years_norm = qdm_matrix_filter(future, FUTURE_MONTH, e->parameters.month, FUTURE_YEAR);
    days = qdm_matrix_filter(future, FUTURE_MONTH, e->parameters.month, FUTURE_DAY);

    gsl_vector_minmax(years, &e->years_min, &e->years_max);

    double years_range = e->years_max - e->years_min;

    status = gsl_vector_add_constant(years_norm, -e->years_min);
    if (status != 0) {
      goto cleanup;
    }

    status = gsl_vector_scale(years_norm, 1 / years_range);
    if (status != 0) {
      goto cleanup;
    }
  }

  /* Extract Values */

  {
    values = qdm_matrix_filter(future, FUTURE_MONTH, e->parameters.month, FUTURE_VALUE);
    values_sorted = qdm_vector_sorted(values);
  }

  /* Optimize Knots */

  {
    gsl_vector *middle = qdm_vector_seq(e->parameters.tau_low, e->parameters.tau_high, 0.005);

    double possible_low = e->parameters.tau_low + 0.1;
    double possible_high = e->parameters.tau_high - 0.1;
    double possible_delta = (possible_high - possible_low) / 19;

    gsl_vector *possible_knots = qdm_vector_seq(possible_low, possible_high, possible_delta);

    interior_knots = gsl_vector_calloc(e->parameters.knot);

    status = qdm_knots_optimize(
        interior_knots,
        values_sorted,
        middle,
        possible_knots,

        e->parameters.knot_try,
        e->parameters.spline_df
    );
    if (status != 0) {
      goto cleanup;
    }

    gsl_vector_free(possible_knots);
    gsl_vector_free(middle);
  }

  /* Optimize Theta */

  {
    size_t m = e->parameters.spline_df + e->parameters.knot;

    m_knots = qdm_knots_vector(e->parameters.spline_df, interior_knots);
    gsl_vector *middle = qdm_vector_seq(e->parameters.tau_low, e->parameters.tau_high, 0.01);

    gsl_matrix *ix = gsl_matrix_alloc(middle->size, m+1);
    qdm_ispline_matrix(ix, middle, e->parameters.spline_df, m_knots);

    theta = gsl_matrix_alloc(2, ix->size2);
    gsl_vector_view theta0 = gsl_matrix_row(theta, 0);
    gsl_vector_view theta1 = gsl_matrix_row(theta, 1);

    size_t theta_pivot = qdm_vector_greater_than(years_norm, 0.3);
    gsl_vector_view values0 = gsl_vector_subvector(values, 0, theta_pivot);
    gsl_vector_view values1 = gsl_vector_subvector(values, theta_pivot, values->size - theta_pivot);
    assert(values0.vector.size + values1.vector.size == values->size);

    gsl_vector *eq0 = qdm_vector_quantile(&values0.vector, middle);
    status = qdm_theta_optimize(&theta0.vector, eq0, ix);
    if (status != 0) {
      goto cleanup;
    }

    gsl_vector *eq1 = qdm_vector_quantile(&values1.vector, middle);
    status = qdm_theta_optimize(&theta1.vector, eq1, ix);
    if (status != 0) {
      goto cleanup;
    }

    status = gsl_vector_sub(&theta1.vector, &theta0.vector);
    if (status != 0) {
      goto cleanup;
    }

    qdm_theta_matrix_constrain(theta, e->parameters.theta_min);

    ll = gsl_vector_alloc(values->size);
    tau = gsl_vector_alloc(values->size);
   
    for (size_t i = 0; i < values->size; i++) {
      status = qdm_logl_2(
          &ll->data[i],
          &tau->data[i],

          gsl_vector_get(years_norm, i),
          gsl_vector_get(values, i),

          e->parameters.tau_low,
          e->parameters.tau_high,

          e->parameters.xi_low,
          e->parameters.xi_high,

          e->parameters.spline_df,

          theta,
          m_knots
      );
      if (status != 0) {
        goto cleanup;
      }
    }

    gsl_vector_free(eq1);
    gsl_vector_free(eq0);
    gsl_matrix_free(ix);
    gsl_vector_free(middle);
  }

  /* MCMC */

  {
    theta_tune = gsl_matrix_alloc(theta->size1, theta->size2);
    gsl_matrix_set_all(theta_tune, e->parameters.theta_tune_sd);

    theta_acc = gsl_matrix_calloc(theta->size1, theta->size2);

    gsl_matrix *theta_p = gsl_matrix_alloc(theta->size1, theta->size2);
    
    gsl_matrix *theta_star = gsl_matrix_alloc(theta->size1, theta->size2);
    gsl_matrix_memcpy(theta_star, theta);

    gsl_matrix *theta_star_p = gsl_matrix_alloc(theta->size1, theta->size2);

    gsl_vector *xi = gsl_vector_alloc(2);
    gsl_vector_set(xi, 0, e->parameters.xi_low);
    gsl_vector_set(xi, 1, e->parameters.xi_high);

    xi_tune = gsl_vector_alloc(2);
    gsl_vector_set_all(xi_tune, e->parameters.xi_tune_sd);

    xi_acc = gsl_vector_alloc(2);
    gsl_vector_set_all(xi_tune, e->parameters.xi_tune_sd);

    gsl_vector *ll_p = gsl_vector_alloc(values->size);
    gsl_vector *tau_p = gsl_vector_alloc(values->size);

    for (int iteration = 0; iteration < e->parameters.iter + e->parameters.burn; iteration++) {
      if (iteration % 100 == 0) {
        fprintf(stderr, "iteration: %d\n", iteration);
      }

      /* Update the mth spline. */
      for (size_t p = 0; p < theta->size1; p++) {
        for (size_t m = 0; m < theta->size2; m++) {
          /* Propose new theta star. */
          status = gsl_matrix_memcpy(theta_star_p, theta_star);
          if (status != 0) {
            goto cleanup;
          }

          double proposal = gsl_matrix_get(theta_tune, p, m) * gsl_ran_ugaussian(qdm_state->rng) + gsl_matrix_get(theta_star, p, m);
          gsl_matrix_set(theta_star_p, p, m, proposal);

          gsl_matrix_memcpy(theta_p, theta_star_p);
          qdm_theta_matrix_constrain(theta_p, e->parameters.theta_min);

          for (size_t i = 0; i < values->size; i++) {
            status = qdm_logl_2(
                &ll_p->data[i],
                &tau_p->data[i],

                gsl_vector_get(years_norm, i),
                gsl_vector_get(values, i),

                e->parameters.tau_low,
                e->parameters.tau_high,

                gsl_vector_get(xi, 0),
                gsl_vector_get(xi, 1),

                e->parameters.spline_df,

                theta_p,
                m_knots
            );
            if (status != 0) {
              goto cleanup;
            }
          }

          double ratio = qdm_vector_sum(ll_p) - qdm_vector_sum(ll);
          if (log(gsl_rng_uniform(qdm_state->rng)) < ratio) {
            gsl_matrix_set(theta_acc, p, m, gsl_matrix_get(theta_acc, p, m) + 1);

            status = gsl_matrix_memcpy(theta_star, theta_star_p);
            if (status != 0) {
              goto cleanup;
            }

            status = gsl_matrix_memcpy(theta, theta_p);
            if (status != 0) {
              goto cleanup;
            }

            status = gsl_vector_memcpy(ll, ll_p);
            if (status != 0) {
              goto cleanup;
            }
          }
        }
      }

      /* Update xi_low. */
      {
        double xi_low = gsl_vector_get(xi, 0);
        double xi_high = gsl_vector_get(xi, 1);

        double xi_low_p = exp(log(xi_low) + gsl_vector_get(xi_tune, 0) * gsl_ran_ugaussian(qdm_state->rng));

        for (size_t i = 0; i < values->size; i++) {
          status = qdm_logl_2(
              &ll_p->data[i],
              &tau_p->data[i],

              gsl_vector_get(years_norm, i),
              gsl_vector_get(values, i),

              e->parameters.tau_low,
              e->parameters.tau_high,

              xi_low_p,
              xi_high,

              e->parameters.spline_df,

              theta,
              m_knots
          );
          if (status != 0) {
            goto cleanup;
          }
        }

        double ratio = -0.5 * (1 / e->parameters.xi_prior_var) *  pow(log(xi_low_p) - e->parameters.xi_prior_mean, 2) +
                        0.5 * (1 / e->parameters.xi_prior_var) * (pow(log(xi_low  ) - e->parameters.xi_prior_mean, 2) + qdm_vector_sum(ll_p) - qdm_vector_sum(ll));
        if (log(gsl_rng_uniform(qdm_state->rng)) < ratio) {
          gsl_vector_set(xi, 0, xi_low_p);

          status = gsl_vector_memcpy(ll, ll_p);
          if (status != 0) {
            goto cleanup;
          }

          gsl_vector_set(xi_acc, 0, gsl_vector_get(xi_acc, 0) + 1);
        }
      }

      /* Update xi_high. */
      {
        double xi_low = gsl_vector_get(xi, 0);
        double xi_high = gsl_vector_get(xi, 1);

        double xi_high_p = exp(log(xi_high) + gsl_vector_get(xi_tune, 1) * gsl_ran_ugaussian(qdm_state->rng));
        for (size_t i = 0; i < values->size; i++) {
          status = qdm_logl_2(
              &ll_p->data[i],
              &tau_p->data[i],

              gsl_vector_get(years_norm, i),
              gsl_vector_get(values, i),

              e->parameters.tau_low,
              e->parameters.tau_high,

              xi_low,
              xi_high_p,

              e->parameters.spline_df,

              theta,
              m_knots
          );
          if (status != 0) {
            goto cleanup;
          }
        }

        double ratio = -0.5 * (1 / e->parameters.xi_prior_var) *  pow(log(xi_high_p) - e->parameters.xi_prior_mean, 2) +
                        0.5 * (1 / e->parameters.xi_prior_var) * (pow(log(xi_high  ) - e->parameters.xi_prior_mean, 2) + qdm_vector_sum(ll_p) - qdm_vector_sum(ll));
        if (log(gsl_rng_uniform(qdm_state->rng)) < ratio) {
          gsl_vector_set(xi, 1, xi_high_p);

          status = gsl_vector_memcpy(ll, ll_p);
          if (status != 0) {
            goto cleanup;
          }

          gsl_vector_set(xi_acc, 1, gsl_vector_get(xi_acc, 1) + 1);
        }
      }

      /* Keep the keepers. */
      if (iteration > e->parameters.burn && iteration % e->parameters.thin == 0) {
        size_t idx = (iteration - e->parameters.burn) / e->parameters.thin - 1;
        e->intermediates[idx] = (qdm_intermediate_result) {
          .theta      = qdm_matrix_copy(theta),
          .theta_star = qdm_matrix_copy(theta_star),
          .ll         = qdm_vector_copy(ll),
          .tau        = qdm_vector_copy(tau),
          .xi         = qdm_vector_copy(xi),
        };
      }

      /* Update tuning parameters. */
      if (iteration < e->parameters.burn && (iteration + 1) % e->parameters.acc_check == 0) {
        for (size_t p = 0; p < theta->size1; p++) {
          for (size_t m = 0; m < theta->size2; m++) {
            if (gsl_matrix_get(theta_acc, p, m) / (double)e->parameters.acc_check > 0.5) {
              gsl_matrix_set(theta_tune, p, m, gsl_min(gsl_matrix_get(theta_tune, p, m) * 1.2, 5));
            }

            if (gsl_matrix_get(theta_acc, p, m) / (double)e->parameters.acc_check < 0.3) {
              gsl_matrix_set(theta_tune, p, m, gsl_matrix_get(theta_tune, p, m) * 0.8);
            }
          }
        }

        if (gsl_vector_get(xi_acc, 0) / e->parameters.acc_check > 0.5) {
          gsl_vector_set(xi_tune, 0, gsl_min(gsl_vector_get(xi_tune, 0) * 1.2, 5));
        }
        if (gsl_vector_get(xi_acc, 0) / e->parameters.acc_check < 0.3) {
          gsl_vector_set(xi_tune, 0, gsl_vector_get(xi_tune, 0) * 0.8);
        }

        if (gsl_vector_get(xi_acc, 1) / e->parameters.acc_check > 0.5) {
          gsl_vector_set(xi_tune, 1, gsl_min(gsl_vector_get(xi_tune, 1) * 1.2, 5));
        }
        if (gsl_vector_get(xi_acc, 1) / e->parameters.acc_check < 0.3) {
          gsl_vector_set(xi_tune, 1, gsl_vector_get(xi_tune, 1) * 0.8);
        }

        /* Reset the acceptance counter. */
        gsl_matrix_set_zero(theta_acc);
        gsl_vector_set_zero(xi_acc);
      }
    }

    gsl_vector_free(tau_p);
    gsl_vector_free(ll_p);
    gsl_vector_free(xi);
    gsl_matrix_free(theta_star_p);
    gsl_matrix_free(theta_star);
    gsl_matrix_free(theta_p);
  }

  /* Model Diagnostics */

  {
    // WAIC

    gsl_vector *ll_tmp = gsl_vector_alloc(e->intermediates_count);
    gsl_vector *ll_exp = gsl_vector_alloc(e->intermediates_count);

    gsl_vector *lppd = gsl_vector_alloc(values->size);
    gsl_vector *vlog = gsl_vector_alloc(values->size);

    for (size_t i = 0; i < values->size; i++) {
      for (size_t j = 0; j < e->intermediates_count; j++) {
        double tau_tmp = 0;

        status = qdm_logl_2(
            &ll_tmp->data[j],
            &tau_tmp,

            gsl_vector_get(years_norm, i),
            gsl_vector_get(values, i),

            e->parameters.tau_low,
            e->parameters.tau_high,

            gsl_vector_get(e->intermediates[j].xi, 0),
            gsl_vector_get(e->intermediates[j].xi, 1),

            e->parameters.spline_df,

            e->intermediates[j].theta,
            m_knots
        );
        if (status != 0) {
          goto cleanup;
        }

        gsl_vector_set(ll_exp, j, exp(gsl_vector_get(ll_tmp, j)));
      }

      gsl_vector_set(lppd, i, log(gsl_stats_mean(ll_exp->data, ll_exp->stride, ll_exp->size)));
      gsl_vector_set(vlog, i, gsl_stats_variance(ll_tmp->data, ll_tmp->stride, ll_tmp->size));
    }

    status = gsl_vector_sub(lppd, vlog);
    if (status != 0) {
      goto cleanup;
    }

    e->waic = -2 * qdm_vector_sum(lppd);
    e->pwaic = qdm_vector_sum(vlog);

    // DIC

    gsl_vector *ll_sum = gsl_vector_alloc(e->intermediates_count);

    for (size_t i = 0; i < e->intermediates_count; i++) {
      gsl_vector_set(ll_sum, i, qdm_vector_sum(e->intermediates[i].ll));
    }

    status = gsl_vector_scale(ll_sum, -2);
    if (status != 0) {
      goto cleanup;
    }

    double d_bar = gsl_stats_mean(ll_sum->data, ll_sum->stride, ll_sum->size);
    double d_theta_bar = 0;

    gsl_vector *xi_low = gsl_vector_alloc(e->intermediates_count);
    gsl_vector *xi_high = gsl_vector_alloc(e->intermediates_count);

    for (size_t i = 0; i < e->intermediates_count; i++) {
      gsl_vector_set(xi_low, i, gsl_vector_get(e->intermediates[i].xi, 0));
      gsl_vector_set(xi_high, i, gsl_vector_get(e->intermediates[i].xi, 1));
    }

    double xi_low_bar = gsl_stats_mean(xi_low->data, xi_low->stride, xi_low->size);
    double xi_high_bar = gsl_stats_mean(xi_high->data, xi_high->stride, xi_high->size);

    gsl_vector *theta_ij = gsl_vector_alloc(e->intermediates_count);
    gsl_matrix *theta_bar = gsl_matrix_alloc(theta->size1, theta->size2);

    for (size_t i = 0; i < theta->size1; i++) {
      for (size_t j = 0; j < theta->size2; j++) {
        for (size_t k = 0; k < e->intermediates_count; k++) {
          gsl_vector_set(theta_ij, k, gsl_matrix_get(e->intermediates[k].theta, i, j));
        }

        gsl_matrix_set(theta_bar, i, j, gsl_stats_mean(theta_ij->data, theta_ij->stride, theta_ij->size));
      }
    }

    for (size_t i = 0; i < values->size; i++) {
      double theta_bar_ll = 0;
      double tau_tmp = 0;

      status = qdm_logl_2(
          &theta_bar_ll,
          &tau_tmp,

          gsl_vector_get(years_norm, i),
          gsl_vector_get(values, i),

          e->parameters.tau_low,
          e->parameters.tau_high,

          xi_low_bar,
          xi_high_bar,

          e->parameters.spline_df,

          theta_bar,
          m_knots
      );
      if (status != 0) {
        goto cleanup;
      }

      d_theta_bar = d_theta_bar + -2 * theta_bar_ll;
    }

    e->pd = d_bar - d_theta_bar;
    e->dic = e->pd + d_bar;

    e->theta_bar = qdm_matrix_copy(theta_bar);
    e->xi_low_bar = xi_low_bar;
    e->xi_high_bar = xi_high_bar;

    gsl_matrix_free(theta_bar);
    gsl_vector_free(theta_ij);
    gsl_vector_free(xi_high);
    gsl_vector_free(xi_low);
    gsl_vector_free(ll_sum);

    gsl_vector_free(vlog);
    gsl_vector_free(lppd);
    gsl_vector_free(ll_exp);
    gsl_vector_free(ll_tmp);
  }

  /* Bias Correction */

  {
    gsl_vector *h_os = NULL;
    gsl_vector *h_os_sorted = NULL;

    gsl_vector *h_vs = NULL;
    gsl_vector *h_vs_sorted = NULL;

    gsl_vector *bins = NULL;

    gsl_vector *h_os_qs = NULL;
    gsl_vector *h_vs_qs = NULL;

    gsl_vector *track = NULL;

    h_os = qdm_matrix_filter(
        historical,
        HISTORICAL_MONTH,
        e->parameters.month,
        HISTORICAL_OBSERVED
    );
    h_os_sorted = qdm_vector_sorted(h_os);

    h_vs = qdm_matrix_filter(
        historical,
        HISTORICAL_MONTH,
        e->parameters.month,
        HISTORICAL_VALUE
    );
    h_vs_sorted = qdm_vector_sorted(h_vs);

    bins = qdm_vector_seq(0, 1, 0.0001);

    h_os_qs = qdm_vector_quantile(h_os_sorted, bins);
    h_vs_qs = qdm_vector_quantile(h_vs_sorted, bins);

    track = gsl_vector_alloc(e->intermediates_count);

    e->bias_corrected = gsl_matrix_alloc(values->size, 5);

    for (size_t i = 0; i < values->size; i++) {
      for (size_t j = 0; j < e->intermediates_count; j++) {
        double ll_tmp = 0;
        double tau_tmp = 0;

        status = qdm_logl_2(
            &ll_tmp,
            &tau_tmp,

            gsl_vector_get(years_norm, i),
            gsl_vector_get(values, i),

            e->parameters.tau_low,
            e->parameters.tau_high,

            gsl_vector_get(e->intermediates[j].xi, 0),
            gsl_vector_get(e->intermediates[j].xi, 1),

            e->parameters.spline_df,

            e->intermediates[j].theta,
            m_knots
        );
        if (status != 0) {
          goto cleanup;
        }

        size_t tau_idx = 0;
        for (; tau_idx < bins->size; tau_idx++) {
          if (tau_tmp < gsl_vector_get(bins, tau_idx)) {
            break;
          }
        }

        double c1 = gsl_vector_get(h_os_qs, tau_idx);
        double c2 = gsl_vector_get(values, i) - gsl_vector_get(h_vs_qs, tau_idx);

        gsl_vector_set(track, j, c1 + c2);
      }

      double mean = gsl_stats_mean(track->data, track->stride, track->size);
      double sd = gsl_stats_sd(track->data, track->stride, track->size);

      gsl_matrix_set(e->bias_corrected, i, 0, gsl_vector_get(years, i));
      gsl_matrix_set(e->bias_corrected, i, 1, e->parameters.month);
      gsl_matrix_set(e->bias_corrected, i, 2, gsl_vector_get(days, i));
      gsl_matrix_set(e->bias_corrected, i, 3, mean);
      gsl_matrix_set(e->bias_corrected, i, 4, sd);
    }

    gsl_vector_free(h_os);
    gsl_vector_free(h_os_sorted);

    gsl_vector_free(h_vs);
    gsl_vector_free(h_vs_sorted);

    gsl_vector_free(bins);

    gsl_vector_free(h_os_qs);
    gsl_vector_free(h_vs_qs);

    gsl_vector_free(track);
  }

  /* Results */
  {
    e->theta_acc = qdm_matrix_copy(theta_acc);
    e->xi_acc = qdm_vector_copy(xi_acc);

    e->m_knots = qdm_vector_copy(m_knots);
  }

  gsl_vector_free(years);
  gsl_vector_free(years_norm);
  gsl_vector_free(days);
  gsl_vector_free(values);
  gsl_vector_free(values_sorted);
  gsl_vector_free(interior_knots);
  gsl_vector_free(m_knots);
  gsl_matrix_free(theta);
  gsl_vector_free(ll);
  gsl_vector_free(tau);

  gsl_matrix_free(theta_tune);
  gsl_matrix_free(theta_acc);
  gsl_vector_free(xi_tune);
  gsl_vector_free(xi_acc);

cleanup:
  e->elapsed = ((double)(clock() - started_at)) / CLOCKS_PER_SEC;

  return status;
}

int
main(int argc, char **argv)
{
  int status = 0;

  const char *progname = "qdm";
  int nerrors = 0;

  struct arg_file *input_path = NULL;
  struct arg_file *output_path = NULL;
  struct arg_lit *help = NULL;
  struct arg_end *end = NULL;

  void *argtable[] = {
    input_path  = arg_file0("i", "input",  "PATH", "path to input file"),
    output_path = arg_file0("o", "output", "PATH", "path to output file"),

    help        = arg_lit0(NULL, "help", "display this help and exit"),

    end         = arg_end(MAX_ARG_ERRORS),
  };

  hid_t input = -1;

  gsl_matrix *historical = NULL;
  gsl_matrix *future = NULL;

  hid_t output = -1;

  if (arg_nullcheck(argtable) != 0) {
    fprintf(stderr, "%s: insufficient memory\n", progname);

    status = 1;
    goto cleanup;
  }

  /* Set default values. */

  input_path->filename[0]  = "input.hd5";
  output_path->filename[0] = "output.hd5";

  /* Parse command line. */

  nerrors = arg_parse(argc, argv, argtable);

  if (help->count > 0) {
    printf("Usage: %s", progname);
    arg_print_syntax(stdout, argtable, "\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");

    status = 0;
    goto cleanup;
  }

  if (nerrors > 0) {
    arg_print_errors(stderr, end, progname);
    fprintf(stderr, "Try '%s --help' for more information.\n", progname);

    status = 1;
    goto cleanup;
  }

  /* Load Config */

  qdm_config *cfg = NULL;
  status = qdm_config_read(&cfg, input_path->filename[0], "input/config");
  if (status != 0) {
    fprintf(stderr, "Failed to read config from input file.\n");

    goto cleanup;
  }

  fprintf(stderr, "config:\n");
  qdm_config_fwrite(stderr, cfg);

  /* Load Input Data */

  input = H5Fopen(input_path->filename[0], H5F_ACC_RDONLY, H5P_DEFAULT);
  if (input < 0) {
    fprintf(stderr, "Failed to open input file.\n");

    status = input;
    goto cleanup;
  }

  status = qdm_matrix_hd5_read(input, "input/data/historical", &historical);
  if (status != 0) {
    fprintf(stderr, "Failed to read historical input matrix.\n");

    goto cleanup;
  }

  status = qdm_matrix_hd5_read(input, "input/data/future", &future);
  if (status != 0) {
    fprintf(stderr, "Failed to read future input matrix.\n");

    goto cleanup;
  }

  /* Setup Output */

  output = qdm_data_create_file(output_path->filename[0]);
  if (output < 0) {
    status = -1;
    goto cleanup;
  }

  /* Perform Evaluations */

  size_t evaluations_expected = cfg->months->size * cfg->knots->size * cfg->tau_highs->size * cfg->tau_lows->size;
  size_t evaluations_done = 0;
  double evaluations_elapsed = 0;

  for (size_t months_idx = 0; months_idx < cfg->months->size; months_idx++) {
    for (size_t knots_idx = 0; knots_idx < cfg->knots->size; knots_idx++) {
      for (size_t tau_highs_idx = 0; tau_highs_idx < cfg->tau_highs->size; tau_highs_idx++) {
        for (size_t tau_lows_idx = 0; tau_lows_idx < cfg->tau_lows->size; tau_lows_idx++) {
          qdm_parameters parameters = {
            .acc_check = cfg->acc_check,
            .burn = cfg->burn,
            .iter = cfg->iter,
            .knot_try = cfg->knot_try,
            .spline_df = cfg->spline_df,
            .thin = cfg->thin,

            .month = gsl_vector_get(cfg->months, months_idx),
            .knot = gsl_vector_get(cfg->knots, knots_idx),
            .tau_high = gsl_vector_get(cfg->tau_highs, tau_highs_idx),
            .tau_low = gsl_vector_get(cfg->tau_lows, tau_lows_idx),

            .theta_min = cfg->theta_min,
            .theta_tune_sd = cfg->theta_tune_sd,

            .xi_high = cfg->xi_high,
            .xi_low = cfg->xi_low,
            .xi_prior_mean = cfg->xi_prior_mean,
            .xi_prior_var = cfg->xi_prior_var,
            .xi_tune_sd = cfg->xi_tune_sd,
          };

          fprintf(stderr, "parameters:\n");
          qdm_parameters_fprint(stderr, &parameters);

          qdm_evaluation *evaluation = qdm_evaluation_new(&parameters);

          status = qdm_evaluation_run(evaluation, historical, future);
          if (status != 0) {
            goto cleanup;
          }

          evaluations_done += 1;
          evaluations_elapsed += evaluation->elapsed;

          fprintf(stderr, "evaluation[%zu/%zu]:\n", evaluations_done, evaluations_expected);
          qdm_evaluation_fprint(stderr, evaluation);

          double spe = evaluations_elapsed / (double)evaluations_done;
          double rem = (double)(evaluations_expected - evaluations_done);
          double eta = rem * spe;
          fprintf(stderr, "eta: %f\n", eta);

          {
            char group_name[PATH_MAX];

            snprintf(group_name, PATH_MAX, "/output/evaluation/e%0*zu", (int)(evaluations_expected % 10), evaluations_done);

            hid_t group = qdm_data_create_group(output, group_name);
            if (group < 0) {
              status = group;
              goto cleanup;
            }

            status = qdm_evaluation_write(group, evaluation);
            if (status != 0) {
              goto cleanup;
            }

            H5Gclose(group);

            status = H5Fflush(output, H5F_SCOPE_GLOBAL);
            if (status != 0) {
              goto cleanup;
            }
          }

          qdm_evaluation_free(evaluation);

          sleep(5);
        }
      }
    }
  }

cleanup:
  if (output >= 0) {
    H5Fclose(output);
  }

  gsl_matrix_free(future);
  gsl_matrix_free(historical);

  if (input >= 0) {
    H5Fclose(input);
  }

  arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

  return -status;
}
