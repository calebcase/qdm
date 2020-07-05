#include "qdm.h"

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

  WRITE_DOUBLE(bound);

#undef WRITE_DOUBLE
#undef WRITE_INT

  return status;
}
