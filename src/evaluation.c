#include <assert.h>

#include <gsl/gsl_statistics_double.h>

#include <qdm.h>

qdm_evaluation *
qdm_evaluation_new(const qdm_parameters *parameters)
{
  qdm_evaluation *e = malloc(sizeof(qdm_evaluation));

  e->parameters = *parameters;

  e->years_min = 0;
  e->years_max = 0;

  e->waic = 0;
  e->pwaic = 0;
  e->dic = 0;
  e->pd = 0;

  e->mcmc = NULL;

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
  if (e == NULL) {
    return;
  }

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

  qdm_mcmc_free(e->mcmc);
  e->mcmc = NULL;

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

  fprintf(f, "%swaic: %f\n", prefix, e->waic);
  fprintf(f, "%spwaic: %f\n", prefix, e->pwaic);
  fprintf(f, "%sdic: %f\n", prefix, e->dic);
  fprintf(f, "%spd: %f\n", prefix, e->pd);

  fprintf(f, "%smcmc r.s: %zu\n", prefix, e->mcmc->r.s);

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
  hid_t mcmc_group = -1;

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

  WRITE_DOUBLE(waic);
  WRITE_DOUBLE(pwaic);
  WRITE_DOUBLE(dic);
  WRITE_DOUBLE(pd);

  mcmc_group = qdm_data_create_group(id, "mcmc");
  if (mcmc_group < 0) {
    status = mcmc_group;

    goto cleanup;
  }

  status = qdm_mcmc_write(mcmc_group, e->mcmc);
  if (status != 0) {
    goto cleanup;
  }

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

  if (mcmc_group >= 0) {
    H5Gclose(mcmc_group);
  }

  return status;
}

int
qdm_evaluation_run(
    qdm_evaluation *e,
    const gsl_matrix *data,
    size_t data_year_idx,
    size_t data_month_idx,
    size_t data_value_idx
)
{
  int status = 0;

  clock_t started_at = clock();

  gsl_vector *years = NULL;
  gsl_vector *years_norm = NULL;

  gsl_vector *values = NULL;
  gsl_vector *values_sorted = NULL;

  gsl_vector *interior_knots = NULL;

  gsl_vector *m_knots = NULL;
  gsl_matrix *theta = NULL;
  gsl_vector *ll = NULL;
  gsl_vector *tau = NULL;
  qdm_tau *t = NULL;

  gsl_matrix *theta_tune = NULL;
  gsl_matrix *theta_acc = NULL;
  gsl_vector *xi_tune = NULL;
  gsl_vector *xi_acc = NULL;

  gsl_vector *xi = gsl_vector_alloc(2);
  gsl_vector_set(xi, 0, e->parameters.xi_low);
  gsl_vector_set(xi, 1, e->parameters.xi_high);

  /* Initialize RNG */

  gsl_rng_set(qdm_state->rng, e->parameters.rng_seed);

  /* Normalize Years */

  {
    years = qdm_matrix_filter(data, data_month_idx, e->parameters.month, data_year_idx);
    years_norm = qdm_matrix_filter(data, data_month_idx, e->parameters.month, data_year_idx);

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
    values = qdm_matrix_filter(data, data_month_idx, e->parameters.month, data_value_idx);
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

        qdm_state->rng,

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
    t = qdm_tau_alloc(e->parameters.tau_low, e->parameters.tau_high, e->parameters.spline_df, m_knots);

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

    gsl_vector *eq1 = NULL;

    if (e->parameters.truncate) {
      gsl_vector_set_zero(&theta1.vector);
    } else {
      eq1 = qdm_vector_quantile(&values1.vector, middle);
      status = qdm_theta_optimize(&theta1.vector, eq1, ix);
      if (status != 0) {
        goto cleanup;
      }

      status = gsl_vector_sub(&theta1.vector, &theta0.vector);
      if (status != 0) {
        goto cleanup;
      }
    }

    qdm_theta_matrix_constrain(theta, e->parameters.theta_min);

    ll = gsl_vector_alloc(values->size);
    tau = gsl_vector_alloc(values->size);

    for (size_t i = 0; i < values->size; i++) {
      qdm_logl_3(
          &ll->data[i],
          &tau->data[i],

          gsl_vector_get(years_norm, i),
          gsl_vector_get(values, i),

          t,
          xi,

          theta
      );
    }

    gsl_vector_free(eq1);
    gsl_vector_free(eq0);
    gsl_matrix_free(ix);
    gsl_vector_free(middle);
  }

  /* MCMC */

  {
    qdm_mcmc_parameters mcmc_p = {
      .rng = qdm_state->rng,

      .burn = e->parameters.burn,
      .iter = e->parameters.iter,
      .thin = e->parameters.thin,

      .acc_check = e->parameters.acc_check,
      .spline_df = e->parameters.spline_df,

      .xi    = xi,
      .ll    = ll,
      .tau   = tau,
      .theta = theta,

      .xi_prior_mean = e->parameters.xi_prior_mean,
      .xi_prior_var  = e->parameters.xi_prior_var,
      .xi_tune_sd    = e->parameters.xi_tune_sd,

      .theta_min     = e->parameters.theta_min,
      .theta_tune_sd = e->parameters.theta_tune_sd,

      .x = years_norm,
      .y = values,

      .t = t,
      .m_knots = m_knots,

      .truncate = e->parameters.truncate,
    };

    e->mcmc = qdm_mcmc_alloc(mcmc_p);

    status = qdm_mcmc_run(e->mcmc);
    if (status != 0) {
      goto cleanup;
    }
  }

  /* Model Diagnostics */

  {
    // WAIC

    gsl_vector *ll_exp = gsl_vector_calloc(e->mcmc->r.s);

    gsl_vector *lppd = gsl_vector_calloc(values->size);
    gsl_vector *vlog = gsl_vector_calloc(values->size);

    for (size_t j = 0; j < values->size; j++) {
      gsl_vector_view ll_k = qdm_ijk_get_k(e->mcmc->r.ll, 0, j);

      for (size_t i = 0; i < ll_k.vector.size; i++) {
        gsl_vector_set(ll_exp, i, exp(gsl_vector_get(&ll_k.vector, i)));
      }

      gsl_vector_set(lppd, j, log(gsl_stats_mean(ll_exp->data, ll_exp->stride, ll_exp->size)));
      gsl_vector_set(vlog, j, gsl_stats_variance(ll_k.vector.data, ll_k.vector.stride, ll_k.vector.size));
    }

    status = gsl_vector_sub(lppd, vlog);
    if (status != 0) {
      goto cleanup;
    }

    e->waic = -2 * qdm_vector_sum(lppd);
    e->pwaic = qdm_vector_sum(vlog);

    // DIC

    gsl_vector *ll_sum = gsl_vector_alloc(e->mcmc->r.s);

    for (size_t k = 0; k < e->mcmc->r.s; k++) {
      gsl_matrix_view ll_ij_m = qdm_ijk_get_ij(e->mcmc->r.ll, k);
      gsl_vector_view ll_ij = gsl_matrix_row(&ll_ij_m.matrix, 0);

      gsl_vector_set(ll_sum, k, qdm_vector_sum(&ll_ij.vector));
    }

    status = gsl_vector_scale(ll_sum, -2);
    if (status != 0) {
      goto cleanup;
    }

    double d_bar = gsl_stats_mean(ll_sum->data, ll_sum->stride, ll_sum->size);
    double d_theta_bar = 0;

    gsl_vector_view xi_low = qdm_ijk_get_k(e->mcmc->r.xi, 0, 0);
    gsl_vector_view xi_high = qdm_ijk_get_k(e->mcmc->r.xi, 0, 1);

    double xi_low_bar = gsl_stats_mean(xi_low.vector.data, xi_low.vector.stride, xi_low.vector.size);
    double xi_high_bar = gsl_stats_mean(xi_high.vector.data, xi_high.vector.stride, xi_high.vector.size);

    gsl_vector *xi_bar = gsl_vector_alloc(2);
    gsl_vector_set(xi_bar, 0, xi_low_bar);
    gsl_vector_set(xi_bar, 1, xi_high_bar);

    gsl_matrix *theta_bar = gsl_matrix_alloc(theta->size1, theta->size2);

    for (size_t i = 0; i < e->mcmc->r.theta->size1; i++) {
      for (size_t j = 0; j < e->mcmc->r.theta->size2; j++) {
        gsl_vector_view theta_k = qdm_ijk_get_k(e->mcmc->r.theta, i, j);

        gsl_matrix_set(theta_bar, i, j, gsl_stats_mean(theta_k.vector.data, theta_k.vector.stride, theta_k.vector.size));
      }
    }

    for (size_t i = 0; i < values->size; i++) {
      double theta_bar_ll = 0;
      double tau_tmp = 0;

      qdm_logl_3(
          &theta_bar_ll,
          &tau_tmp,

          gsl_vector_get(years_norm, i),
          gsl_vector_get(values, i),

          t,
          xi_bar,

          theta_bar
      );

      d_theta_bar = d_theta_bar + -2 * theta_bar_ll;
    }

    e->pd = d_bar - d_theta_bar;
    e->dic = e->pd + d_bar;

    e->theta_bar = qdm_matrix_copy(theta_bar);
    e->xi_low_bar = xi_low_bar;
    e->xi_high_bar = xi_high_bar;

    gsl_vector_free(xi_bar);
    gsl_matrix_free(theta_bar);
    gsl_vector_free(ll_sum);

    gsl_vector_free(vlog);
    gsl_vector_free(lppd);
    gsl_vector_free(ll_exp);
  }

  /* Results */
  {
    e->theta_acc = qdm_matrix_copy(e->mcmc->w.theta_acc);
    e->xi_acc = qdm_vector_copy(e->mcmc->w.xi_acc);

    e->m_knots = qdm_vector_copy(m_knots);
  }

  gsl_vector_free(years);
  gsl_vector_free(years_norm);
  gsl_vector_free(values);
  gsl_vector_free(values_sorted);
  gsl_vector_free(interior_knots);
  qdm_tau_free(t);
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
