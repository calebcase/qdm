#define _GNU_SOURCE
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <argtable2.h>
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
phase_run(
    qdm_evaluation **result,

    const qdm_config *cfg,
    size_t months_idx,

    const char *phase_name,

    const gsl_matrix *data,
    size_t data_year_idx,
    size_t data_month_idx,
    size_t data_value_idx,

    const bool truncate,

    hid_t output
)
{
  int status = 0;

  /* Attempt to read an existing result. If one isn't found, then proceed with calculating. */
  {
    char group_name[PATH_MAX];

    int month = (int)gsl_vector_get(cfg->months, months_idx);

    snprintf(group_name, PATH_MAX, "/output/m%02d/%s/analysis", month, phase_name);

    hid_t group = H5Gopen(output, group_name, H5P_DEFAULT);
    if (group >= 0) {
      status = qdm_evaluation_read(group, result);
      if (status == 0) {
        fprintf(stderr, "Loaded month %02d phase %s from checkpoint...\n", month, phase_name);

        goto cleanup;
      }
    }

    status = 0;
  }

  qdm_evaluation *best = NULL;

  size_t evaluations_expected = cfg->knots->size * cfg->tau_highs->size * cfg->tau_lows->size;
  size_t evaluations_done = 0;
  double evaluations_elapsed = 0;

  for (size_t knots_idx = 0; knots_idx < cfg->knots->size; knots_idx++) {
    for (size_t tau_highs_idx = 0; tau_highs_idx < cfg->tau_highs->size; tau_highs_idx++) {
      for (size_t tau_lows_idx = 0; tau_lows_idx < cfg->tau_lows->size; tau_lows_idx++) {
        qdm_parameters parameters = {
          .rng_seed = cfg->rng_seed,

          .acc_check = cfg->acc_check,
          .burn = cfg->burn_discovery,
          .iter = cfg->iter_discovery,
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

          .bound = cfg->bound,

          .debug = cfg->debug,
          .tau_table = cfg->tau_table,

          .truncate = truncate,
        };

        fprintf(stderr, "parameters:\n");
        qdm_parameters_fprint(stderr, &parameters);

        qdm_evaluation *evaluation = qdm_evaluation_new(&parameters);

        status = qdm_evaluation_run(evaluation, data, data_year_idx, data_month_idx, data_value_idx);
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

          snprintf(group_name, PATH_MAX, "/output/m%02d/%s/discovery/e%0*zu", (int)parameters.month, phase_name, (int)log10(evaluations_expected), evaluations_done);

          H5Ldelete(output, group_name, H5P_DEFAULT);

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

        if (best == NULL) {
          best = evaluation;
          evaluation = NULL;
        }
        else if (evaluation->waic < best->waic) {
          fprintf(stderr, "Found a better evaluation %f < %f\n", evaluation->waic, best->waic);

          qdm_evaluation_free(best);
          best = evaluation;
          evaluation = NULL;
        }

        qdm_evaluation_free(evaluation);
      }
    }
  }

  fprintf(stderr, "Best Evaluation:\n");
  qdm_evaluation_fprint(stderr, best);

  {
    qdm_parameters parameters = {
      .rng_seed = cfg->rng_seed,

      .acc_check = cfg->acc_check,
      .burn = cfg->burn_analysis,
      .iter = cfg->iter_analysis,
      .knot_try = cfg->knot_try,
      .spline_df = cfg->spline_df,
      .thin = cfg->thin,

      .month = gsl_vector_get(cfg->months, months_idx),
      .knot = best->parameters.knot,

      .tau_high = best->parameters.tau_high,
      .tau_low = best->parameters.tau_low,

      .theta_min = cfg->theta_min,
      .theta_tune_sd = cfg->theta_tune_sd,

      .xi_high = cfg->xi_high,
      .xi_low = cfg->xi_low,
      .xi_prior_mean = cfg->xi_prior_mean,
      .xi_prior_var = cfg->xi_prior_var,
      .xi_tune_sd = cfg->xi_tune_sd,

      .bound = cfg->bound,

      .debug = cfg->debug,
      .tau_table = cfg->tau_table,

      .truncate = truncate,
    };

    fprintf(stderr, "parameters:\n");
    qdm_parameters_fprint(stderr, &parameters);

    qdm_evaluation *evaluation = qdm_evaluation_new(&parameters);

    status = qdm_evaluation_run(evaluation, data, data_year_idx, data_month_idx, data_value_idx);
    if (status != 0) {
      goto cleanup;
    }

    fprintf(stderr, "evaluation final:\n");
    qdm_evaluation_fprint(stderr, evaluation);

    {
      char group_name[PATH_MAX];

      snprintf(group_name, PATH_MAX, "/output/m%02d/%s/analysis", (int)parameters.month, phase_name);

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

    *result = evaluation;
  }

  qdm_evaluation_free(best);

cleanup:
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

  for (size_t months_idx = 0; months_idx < cfg->months->size; months_idx++) {
    int month = (int)gsl_vector_get(cfg->months, months_idx);

    /* Check if the output for this has already been finalized. */
    {
      char group_name[PATH_MAX];
      snprintf(group_name, PATH_MAX, "/output/m%02d", month);

      hid_t group = H5Gopen(output, group_name, H5P_DEFAULT);
      if (group >= 0) {
        int finalized = 0;

        status = H5LTread_dataset_int(group, "finalized", &finalized);
        if (status >= 0) {
          if (finalized) {
            fprintf(stderr, "Skipping %s, output is already finalized.\n", group_name);

            continue;
          }
        }
      }
    }

    /* Phase 1: Future */
    qdm_evaluation *p1e = NULL;
    status = phase_run(
        &p1e,

        cfg,
        months_idx,

        "p1",

        future,
        FUTURE_YEAR,
        FUTURE_MONTH,
        FUTURE_VALUE,

        false,

        output
    );
    if (status != 0) {
      goto cleanup;
    }

    fprintf(stderr, "\n\nPhase 1 Complete\n\n\n");

    /* Phase 2: Historical Observations */
    cfg->rng_seed = p1e->rng_seed;

    qdm_evaluation *p2e = NULL;
    status = phase_run(
        &p2e,

        cfg,
        months_idx,

        "p2",

        historical,
        HISTORICAL_YEAR,
        HISTORICAL_MONTH,
        HISTORICAL_OBSERVED,

        true,

        output
    );
    if (status != 0) {
      goto cleanup;
    }

    fprintf(stderr, "\n\nPhase 2 Complete\n\n\n");

    /* Phase 3: Historical Hindcasts */
    cfg->rng_seed = p2e->rng_seed;

    qdm_evaluation *p3e = NULL;
    status = phase_run(
        &p3e,

        cfg,
        months_idx,

        "p3",

        historical,
        HISTORICAL_YEAR,
        HISTORICAL_MONTH,
        HISTORICAL_VALUE,

        true,

        output
    );
    if (status != 0) {
      goto cleanup;
    }

    fprintf(stderr, "\n\nPhase 3 Complete\n\n\n");

    /* Bias Correction */
    cfg->rng_seed = p3e->rng_seed;

    {
      char group_name[PATH_MAX];
      snprintf(group_name, PATH_MAX, "/output/m%02d/p4", month);

      hid_t group = H5Gopen(output, group_name, H5P_DEFAULT);
      if (group >= 0) {
        fprintf(stderr, "Loaded month %02d phase p4 from checkpoint...\n", month);
      } else {
        gsl_vector *years = qdm_matrix_filter(future, FUTURE_MONTH, (double)month, FUTURE_YEAR);
        gsl_vector *days = qdm_matrix_filter(future, FUTURE_MONTH, (double)month, FUTURE_DAY);
        gsl_vector *y = qdm_matrix_filter(future, FUTURE_MONTH, (double)month, FUTURE_VALUE);

        gsl_matrix *bc = qdm_bias_correct(years, (double)month, days, y, p1e, p2e, p3e);

        group = qdm_data_create_group(output, group_name);

        status = qdm_matrix_hd5_write(group, "bc", bc);
        if (status != 0) {
          goto cleanup;
        }

        H5Gclose(group);

        gsl_matrix_free(bc);
        gsl_vector_free(y);
        gsl_vector_free(days);
        gsl_vector_free(years);
      }
    }

    fprintf(stderr, "\n\nPhase 4 Complete\n\n\n");

    /* Remove debugging and intermediate information. */
    if (!cfg->debug) {
      char group_name[PATH_MAX];

      for (int phase = 1; phase <= 3; phase++) {
        snprintf(group_name, PATH_MAX, "/output/m%02d/p%d/discovery", month, phase);
        H5Ldelete(output, group_name, H5P_DEFAULT);

        snprintf(group_name, PATH_MAX, "/output/m%02d/p%d/analysis/mcmc", month, phase);
        H5Ldelete(output, group_name, H5P_DEFAULT);
      }
    }

    /* Finalize the output file. */
    {
      char path[PATH_MAX];
      snprintf(path, PATH_MAX, "/output/m%02d/finalized", month);

      status = qdm_int_write(output, path, 1);
      if (status != 0) {
        goto cleanup;
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
