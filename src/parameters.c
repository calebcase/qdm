#include <hdf5.h>
#include <hdf5_hl.h>

#include "qdm.h"

void
qdm_parameters_fprint(
    FILE *f,
    const qdm_parameters *p
)
{
  const char *prefix = "    ";

  fprintf(f, "%srng_seed: %lu\n", prefix, p->rng_seed);

  fprintf(f, "%smonth: %f\n", prefix, p->month);
  fprintf(f, "%sknot: %f\n",  prefix, p->knot);

  fprintf(f, "%stau_high: %f\n", prefix, p->tau_high);
  fprintf(f, "%stau_low: %f\n",  prefix, p->tau_low);

  fprintf(f, "%sdebug: %d\n", prefix, p->debug);
  fprintf(f, "%stau_table: %d\n", prefix, p->tau_table);

  fprintf(f, "%struncate: %s\n", prefix, p->truncate ? "true" : "false");
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

  {
    hsize_t dims[1] = {1};
    unsigned long int data[1] = {p->rng_seed};

    status = H5LTmake_dataset(id, "rng_seed", 1, dims, H5T_NATIVE_ULONG, data);
    if (status < 0) {
      goto cleanup;
    }
  }

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

  WRITE_DOUBLE(bound);

#undef WRITE_DOUBLE
#undef WRITE_INT

cleanup:
  return status;
}

int
qdm_parameters_read(
    hid_t id,
    qdm_parameters *p
)
{
  int status = 0;

#define READ(t, n) \
  status = H5LTread_dataset_##t(id, #n, &p->n); \
  if (status < 0) { \
    goto cleanup; \
  }

  status = H5LTread_dataset(id, "rng_seed", H5T_NATIVE_ULONG, &p->rng_seed);
  if (status < 0) {
    goto cleanup;
  }

  READ(int, acc_check);
  READ(int, burn);
  READ(int, iter);
  READ(int, knot_try);
  READ(int, spline_df);
  READ(int, thin);

  READ(double, month);
  READ(double, knot);

  READ(double, tau_high);
  READ(double, tau_low);

  READ(double, theta_min);
  READ(double, theta_tune_sd);

  READ(double, xi_high);
  READ(double, xi_low);
  READ(double, xi_prior_mean);
  READ(double, xi_prior_var);
  READ(double, xi_tune_sd);

  READ(double, bound);

#undef READ

cleanup:
  return status;
}
