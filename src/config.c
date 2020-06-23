/* Functions for loading the config. */

#include <math.h>
#include <stdlib.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <qdm.h>

int
qdm_config_read(
    qdm_config **cfg,
    const char *file_path,
    const char *group_path
)
{
  int status = 0;

  hid_t file_id;
  hid_t group_id;

  // FIXME: Check for negative return indicating failure to open.
  file_id = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);
  group_id = H5Gopen(file_id, group_path, H5P_DEFAULT);

  // FIXME: Check for malloc failure.
  qdm_config *config = malloc(sizeof(qdm_config));

  // Initialize values
#define READ(t, n) \
  status = H5LTread_dataset_##t(group_id, #n, &config->n); \
  if (status < 0) { \
    goto cleanup; \
  }

#define READ_VECTOR(n) \
  status = qdm_vector_hd5_read(group_id, #n, &config->n); \
  if (status < 0) { \
    goto cleanup; \
  }

  READ(int, rng_seed);

  READ(int, acc_check);

  READ(int, burn_discovery);
  READ(int, iter_discovery);

  READ(int, burn_analysis);
  READ(int, iter_analysis);

  READ(int, knot_try);
  READ(int, spline_df);
  READ(int, thin);

  READ_VECTOR(months);
  READ_VECTOR(knots);

  READ_VECTOR(tau_highs);
  READ_VECTOR(tau_lows);

  READ(double, theta_min);
  READ(double, theta_tune_sd);

  READ(double, xi_high);
  READ(double, xi_low);
  READ(double, xi_prior_mean);
  READ(double, xi_prior_var);
  READ(double, xi_tune_sd);

  READ(int, debug);
  READ(int, tau_table);

#undef READ_VECTOR
#undef READ

  *cfg = config;

cleanup:
  if (status < 0) {
    free(config);
    config = NULL;
    *cfg = NULL;
  }

  H5Gclose(group_id);
  H5Fclose(file_id);

  return status;
}

void
qdm_config_fwrite(
    FILE *f,
    const qdm_config *cfg
)
{
  const char *prefix = "  ";

#define PRINT_INT(n) fprintf(f, "%s%s: %d\n", prefix, #n, cfg->n);
#define PRINT_DOUBLE(n) fprintf(f, "%s%s: %f\n", prefix, #n, cfg->n);
#define PRINT_VECTOR(n) \
  fprintf(f, "%s%s:", prefix, #n); \
  \
  if (cfg->n != NULL) { \
    fprintf(f, "\n"); \
    for (size_t i = 0; i < cfg->n->size; i++) { \
      fprintf(f, "%s  - %f\n", prefix, gsl_vector_get(cfg->n, i)); \
    } \
  } \
  else { \
    fprintf(f, " NULL\n"); \
  }

  PRINT_INT(rng_seed);

  PRINT_INT(acc_check);

  PRINT_INT(burn_discovery);
  PRINT_INT(iter_discovery);

  PRINT_INT(burn_analysis);
  PRINT_INT(iter_analysis);

  PRINT_INT(knot_try);
  PRINT_INT(spline_df);
  PRINT_INT(thin);

  PRINT_VECTOR(months);
  PRINT_VECTOR(knots);

  PRINT_VECTOR(tau_highs);
  PRINT_VECTOR(tau_lows);

  PRINT_DOUBLE(theta_min);
  PRINT_DOUBLE(theta_tune_sd);

  PRINT_DOUBLE(xi_high);
  PRINT_DOUBLE(xi_low);
  PRINT_DOUBLE(xi_prior_mean);
  PRINT_DOUBLE(xi_prior_var);
  PRINT_DOUBLE(xi_tune_sd);

  PRINT_INT(debug);
  PRINT_INT(tau_table);

#undef PRINT_VECTOR
#undef PRINT_DOUBLE
#undef PRINT_INT
}
