#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>

#include <argtable2.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_minmax.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_vector.h>
#include <qdm.h>

#include "config.h"

int
qdm_main(qdm_config *cfg)
{
  int status = 0;

  /* Load Data */

  qdm_data_model *records = NULL;
  size_t nrecords = 0;

  qdm_config_fprint(stderr, cfg);

  status = qdm_data_model_read(
      &records,
      &nrecords,
      cfg->data_input,
      "input/model",
      "m1"
  );
  if (status != 0) {
    return status;
  }

  double month = 4;

  gsl_vector *years = qdm_data_model_year_vector(
      records,
      nrecords,
      (qdm_data_model_selector)qdm_data_model_select_by_month,
      &month
  );

  // Normalize years to range 0 to 1.
  double years_min = 0;
  double years_max = 0;
  gsl_vector_minmax(years, &years_min, &years_max);

  double years_range = years_max - years_min;

  status = gsl_vector_add_constant(years, -years_min);
  if (status != 0) {
    return status;
  }

  status = gsl_vector_scale(years, 1 / years_range);
  if (status != 0) {
    return status;
  }

  gsl_vector *values = qdm_data_model_value_vector(
      records,
      nrecords,
      (qdm_data_model_selector)qdm_data_model_select_by_month,
      &month
  );

  for (size_t i = 0; i < values->size; i++) {
    fprintf(stderr, "selected[%zu]: %.17f %.17f\n",
        i,
        gsl_vector_get(years, i),
        gsl_vector_get(values, i)
    );
  }

  /* Optimize Knots */

  gsl_vector *values_sorted = gsl_vector_alloc(values->size);
  status = gsl_vector_memcpy(values_sorted, values);
  if (status != 0) {
    return status;
  }

  gsl_sort_vector(values_sorted);

  gsl_vector *middle = qdm_vector_seq(cfg->tau_low, cfg->tau_high, 0.005);

  double possible_low = cfg->tau_low + 0.1;
  double possible_high = cfg->tau_high - 0.1;
  double possible_delta = (possible_high - possible_low) / 19;

  gsl_vector *possible_knots = qdm_vector_seq(possible_low, possible_high, possible_delta);

  gsl_vector *interior_knots = gsl_vector_calloc(cfg->main_knots);

  status = qdm_knots_optimize(
      interior_knots,
      values_sorted,
      middle,
      possible_knots,

      cfg->main_knot_try,
      cfg->main_spline_df
  );

  for (size_t i = 0; i < interior_knots->size; i++) {
    fprintf(stderr, "interior_knots[%zu]: %.17f\n",
        i,
        gsl_vector_get(interior_knots, i)
    );
  }

  /* Optimize Theta */
  size_t m = cfg->main_spline_df + cfg->main_knots;

  gsl_vector *m_knots = qdm_knots_vector(cfg->main_spline_df, interior_knots);
  middle = qdm_vector_seq(cfg->tau_low, cfg->tau_high, 0.01);

  gsl_matrix *ix = gsl_matrix_alloc(middle->size, m+1);
  qdm_ispline_matrix(ix, middle, cfg->main_spline_df, m_knots);

  gsl_matrix *theta = gsl_matrix_alloc(2, ix->size2);
  gsl_vector_view theta0 = gsl_matrix_row(theta, 0);
  gsl_vector_view theta1 = gsl_matrix_row(theta, 1);

  size_t theta_pivot = qdm_vector_greater_than(years, 0.3);
  gsl_vector_view values0 = gsl_vector_subvector(values, 0, theta_pivot);
  gsl_vector_view values1 = gsl_vector_subvector(values, theta_pivot, values->size - theta_pivot);
  assert(values0.vector.size + values1.vector.size == values->size);

  gsl_vector *eq0 = qdm_vector_quantile(&values0.vector, middle);
  status = qdm_theta_optimize(&theta0.vector, eq0, ix);
  if (status != 0) {
    return status;
  }

  gsl_vector *eq1 = qdm_vector_quantile(&values1.vector, middle);
  status = qdm_theta_optimize(&theta1.vector, eq1, ix);
  if (status != 0) {
    return status;
  }

  status = gsl_vector_sub(&theta1.vector, &theta0.vector);
  if (status != 0) {
    return status;
  }

  qdm_theta_matrix_constrain(theta, cfg->theta_min);
  fprintf(stderr, "theta:\n");
  qdm_matrix_csv_fwrite(stderr, theta);

  gsl_vector *ll = gsl_vector_alloc(values->size);
  gsl_vector *tau = gsl_vector_alloc(values->size);
 
  gsl_vector *mmm = gsl_vector_alloc(theta->size2);
  gsl_vector *x_ith = gsl_vector_alloc(2);
  gsl_vector_set(x_ith, 0, 1);

  for (size_t i = 0; i < values->size; i++) {
    gsl_vector_set(x_ith, 1, gsl_vector_get(years, i));
    fprintf(stderr, "x_ith[%zu]: %.17f %.17f\n",
        i,
        gsl_vector_get(x_ith, 0),
        gsl_vector_get(x_ith, 1)
    );

    status = gsl_blas_dgemv(CblasTrans, 1, theta, x_ith, 0, mmm);
    if (status != 0) {
      return status;
    }
    fprintf(stderr, "mmm:\n");
    qdm_vector_csv_fwrite(stderr, mmm);

    status = qdm_logl(
        &ll->data[i],
        &tau->data[i],

        gsl_vector_get(values, i),
        cfg->main_spline_df,
        m_knots,
        mmm,

        cfg->tau_low,
        cfg->tau_high,
        cfg->xi_low,
        cfg->xi_high
    );
    if (status != 0) {
      return status;
    }
  }

  for (size_t i = 0; i < ll->size; i++) {
    fprintf(stderr, "(ll,tau)[%zu]: %.17f %.17f\n",
        i,
        gsl_vector_get(ll, i),
        gsl_vector_get(tau, i)
    );
  }

  /* MCMC */
  qdm_intermediate_result keep[cfg->main_iter / cfg->main_thin];

  gsl_matrix *theta_tune = gsl_matrix_alloc(theta->size1, theta->size2);
  gsl_matrix_set_all(theta_tune, cfg->theta_tune_sd);

  gsl_matrix *theta_acc = gsl_matrix_calloc(theta->size1, theta->size2);

  gsl_matrix *theta_star = gsl_matrix_alloc(theta->size1, theta->size2);
  gsl_matrix *theta_star_p = gsl_matrix_alloc(theta->size1, theta->size2);

  gsl_vector *xi = gsl_vector_alloc(2);
  gsl_vector_set(xi, 0, cfg->xi_low);
  gsl_vector_set(xi, 1, cfg->xi_high);

  gsl_vector *xi_tune = gsl_vector_alloc(2);
  gsl_vector_set_all(xi_tune, cfg->xi_tune_sd);

  gsl_vector *xi_acc = gsl_vector_alloc(2);
  gsl_vector_set_all(xi_tune, cfg->xi_tune_sd);

  gsl_vector *ll_p = gsl_vector_alloc(values->size);
  gsl_vector *tau_p = gsl_vector_alloc(values->size);

  fprintf(stderr, "main loop: total = %zu\n", cfg->main_iter + cfg->main_burn);
  for (size_t iteration = 0; iteration < cfg->main_iter + cfg->main_burn; iteration++) {
    fprintf(stderr, "main loop: iteration = %zu\n", iteration);

    /* Update the mth spline. */
    for (size_t p = 0; p < theta->size1; p++) {
      for (size_t m = 0; m < theta->size2; m++) {
        /* Propose new theta star. */
        status = gsl_matrix_memcpy(theta_star_p, theta_star);
        if (status != 0) {
          return status;
        }

        double proposal = gsl_matrix_get(theta_tune, p, m) * gsl_ran_ugaussian(qdm_state->rng) + gsl_matrix_get(theta_star, p, m);
        gsl_matrix_set(theta_star_p, p, m, proposal);

        qdm_theta_matrix_constrain(theta_star_p, cfg->theta_min);

        for (size_t i = 0; i < values->size; i++) {
          gsl_vector_set(x_ith, 1, gsl_vector_get(years, i));

          status = gsl_blas_dgemv(CblasTrans, 1, theta_star_p, x_ith, 0, mmm);
          if (status != 0) {
            return status;
          }

          status = qdm_logl(
              &ll_p->data[i],
              &tau_p->data[i],

              gsl_vector_get(values, i),
              cfg->main_spline_df,
              m_knots,
              mmm,

              cfg->tau_low,
              cfg->tau_high,
              gsl_vector_get(xi, 0),
              gsl_vector_get(xi, 1)
          );
          if (status != 0) {
            return status;
          }
        }

        double ratio = qdm_vector_sum(ll_p) - qdm_vector_sum(ll);
        if (log(gsl_rng_uniform(qdm_state->rng)) < ratio) {
          gsl_matrix_set(theta_acc, p, m, gsl_matrix_get(theta_acc, p, m) + 1);

          status = gsl_matrix_memcpy(theta_star, theta_star_p);
          if (status != 0) {
            return status;
          }

          status = gsl_matrix_memcpy(theta, theta_star);
          if (status != 0) {
            return status;
          }

          status = gsl_vector_memcpy(ll, ll_p);
          if (status != 0) {
            return status;
          }
        }
      }
    }

    /* Update xi_low. */
    {
      double xi_low = gsl_vector_get(xi, 0);
      double xi_high = gsl_vector_get(xi, 1);

      double xi_low_p = exp(log(xi_low) + gsl_vector_get(xi_tune, 0) * gsl_ran_ugaussian(qdm_state->rng));

      // FIXME: This should be a function...
      for (size_t i = 0; i < values->size; i++) {
        gsl_vector_set(x_ith, 1, gsl_vector_get(years, i));

        status = gsl_blas_dgemv(CblasTrans, 1, theta, x_ith, 0, mmm);
        if (status != 0) {
          return status;
        }

        status = qdm_logl(
            &ll_p->data[i],
            &tau_p->data[i],

            gsl_vector_get(values, i),
            cfg->main_spline_df,
            m_knots,
            mmm,

            cfg->tau_low,
            cfg->tau_high,
            xi_low_p,
            xi_high
        );
        if (status != 0) {
          return status;
        }
      }

      double ratio = -0.5 * (1 / cfg->xi_prior_var) *  pow(log(xi_low_p) - cfg->xi_prior_mean, 2) +
                      0.5 * (1 / cfg->xi_prior_var) * (pow(log(xi_low  ) - cfg->xi_prior_mean, 2) + qdm_vector_sum(ll_p) - qdm_vector_sum(ll));
      if (log(gsl_rng_uniform(qdm_state->rng)) < ratio) {
        gsl_vector_set(xi, 0, xi_low_p);

        status = gsl_vector_memcpy(ll, ll_p);
        if (status != 0) {
          return status;
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
        gsl_vector_set(x_ith, 1, gsl_vector_get(years, i));

        status = gsl_blas_dgemv(CblasTrans, 1, theta, x_ith, 0, mmm);
        if (status != 0) {
          return status;
        }

        status = qdm_logl(
            &ll_p->data[i],
            &tau_p->data[i],

            gsl_vector_get(values, i),
            cfg->main_spline_df,
            m_knots,
            mmm,

            cfg->tau_low,
            cfg->tau_high,
            xi_low,
            xi_high_p
        );
        if (status != 0) {
          return status;
        }
      }

      double ratio = -0.5 * (1 / cfg->xi_prior_var) *  pow(log(xi_high_p) - cfg->xi_prior_mean, 2) +
                      0.5 * (1 / cfg->xi_prior_var) * (pow(log(xi_high  ) - cfg->xi_prior_mean, 2) + qdm_vector_sum(ll_p) - qdm_vector_sum(ll));
      if (log(gsl_rng_uniform(qdm_state->rng)) < ratio) {
        gsl_vector_set(xi, 1, xi_high_p);

        status = gsl_vector_memcpy(ll, ll_p);
        if (status != 0) {
          return status;
        }

        gsl_vector_set(xi_acc, 1, gsl_vector_get(xi_acc, 1) + 1);
      }
    }

    /* Keep the keepers. */
    if (iteration > cfg->main_burn && iteration % cfg->main_thin == 0) {
      size_t idx = (iteration - cfg->main_burn) / cfg->main_thin;
      keep[idx] = (qdm_intermediate_result) {
        .theta      = qdm_matrix_copy(theta),
        .theta_star = qdm_matrix_copy(theta_star),
        .ll         = qdm_vector_copy(ll),
        .tau        = qdm_vector_copy(tau),
        .xi         = qdm_vector_copy(xi),
      };

      char *group_path = NULL;
      status = asprintf(&group_path, "output/intermediate/iteration_%zu", iteration);
      if (status < 0) {
        return status;
      }

      status = qdm_data_intermediate_result_write(cfg->data_output, group_path, &keep[idx]);
      if (status != 0) {
        return status;
      }
    }

    /* Update tuning parameters. */
    if (iteration < cfg->main_burn && iteration % cfg->main_acc_check == 0) {
      for (size_t p = 0; p < theta->size1; p++) {
        for (size_t m = 0; m < theta->size2; m++) {
          if (gsl_matrix_get(theta_acc, p, m) / cfg->main_acc_check > 0.5) {
            gsl_matrix_set(theta_tune, p, m, gsl_matrix_get(theta_tune, p, m) * 1.2);
          }

          if (gsl_matrix_get(theta_acc, p, m) / cfg->main_acc_check < 0.3) {
            gsl_matrix_set(theta_tune, p, m, gsl_matrix_get(theta_tune, p, m) * 0.8);
          }
        }
      }

      if (gsl_vector_get(xi_acc, 0) / cfg->main_acc_check > 0.5) {
        gsl_vector_set(xi_tune, 0, gsl_min(gsl_vector_get(xi_tune, 0) * 1.2, 5));
      }
      if (gsl_vector_get(xi_acc, 0) / cfg->main_acc_check < 0.3) {
        gsl_vector_set(xi_tune, 0, gsl_vector_get(xi_tune, 0) * 0.8);
      }

      if (gsl_vector_get(xi_acc, 1) / cfg->main_acc_check > 0.5) {
        gsl_vector_set(xi_tune, 1, gsl_min(gsl_vector_get(xi_tune, 1) * 1.2, 5));
      }
      if (gsl_vector_get(xi_acc, 1) / cfg->main_acc_check < 0.3) {
        gsl_vector_set(xi_tune, 1, gsl_vector_get(xi_tune, 1) * 0.8);
      }

      /* Reset the acceptance counter. */
      gsl_matrix_set_zero(theta_acc);
      gsl_vector_set_zero(xi_acc);
    }
  }

  for (size_t i = 0; i < cfg->main_iter / cfg->main_thin; i++) {
    fprintf(stderr, "keep[%zu]\n", i);
  }

  return status;
}

int
main(int argc, char **argv)
{
  int status = 0;

  const char *progname = "qdm";
  int nerrors = 0;

  struct arg_file *config_path;
  struct arg_lit *help;
  struct arg_end *end;

  void *argtable[] = {
    config_path = arg_file0("c", "config", "PATH", "path to config"),
    help        = arg_lit0(NULL, "help", "display this help and exit"),
    end         = arg_end(20),
  };

  if (arg_nullcheck(argtable) != 0) {
    fprintf(stderr, "%s: insufficient memory\n",progname);

    status = 1;
    goto cleanup;
  }

  /* Set default values. */
  config_path->filename[0] = "qdm.ini";

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

  FILE *config_file = fopen(config_path->filename[0], "r");
  if (config_file == NULL) {
    perror(config_path->filename[0]);

    status = 1;
    goto cleanup;
  }

  qdm_config cfg;
  status = qdm_config_new(config_file, &cfg);
  if (status != 0) {
    fprintf(stderr, "Failed to parse config file.\n");

    status = 1;
    goto cleanup;
  }

  status = qdm_main(&cfg);
  if (status != 0) {
    fprintf(stderr, "Failed to execute qdm routine.\n");

    goto cleanup;
  }

cleanup:
  arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

  return -status;
}
