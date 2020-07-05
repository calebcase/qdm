#include <qdm.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <hdf5.h>
#include <hdf5_hl.h>

static
void
qdm_mcmc_workspace_init(
    qdm_mcmc *mcmc
)
{
  mcmc->w.xi      = gsl_vector_calloc(mcmc->p.xi->size);
  mcmc->w.xi_acc  = gsl_vector_calloc(mcmc->p.xi->size);
  mcmc->w.xi_p    = gsl_vector_calloc(mcmc->p.xi->size);
  mcmc->w.xi_tune = gsl_vector_calloc(mcmc->p.xi->size);

  gsl_vector_memcpy(mcmc->w.xi,   mcmc->p.xi);
  gsl_vector_memcpy(mcmc->w.xi_p, mcmc->p.xi);

  gsl_vector_set_all(mcmc->w.xi_tune, mcmc->p.xi_tune_sd);

  mcmc->w.ll   = gsl_vector_calloc(mcmc->p.ll->size);
  mcmc->w.ll_p = gsl_vector_calloc(mcmc->p.ll->size);

  gsl_vector_memcpy(mcmc->w.ll,   mcmc->p.ll);
  gsl_vector_memcpy(mcmc->w.ll_p, mcmc->p.ll);

  mcmc->w.tau   = gsl_vector_calloc(mcmc->p.tau->size);
  mcmc->w.tau_p = gsl_vector_calloc(mcmc->p.tau->size);

  gsl_vector_memcpy(mcmc->w.tau,   mcmc->p.tau);
  gsl_vector_memcpy(mcmc->w.tau_p, mcmc->p.tau);

  mcmc->w.theta        = gsl_matrix_calloc(mcmc->p.theta->size1, mcmc->p.theta->size2);
  mcmc->w.theta_acc    = gsl_matrix_calloc(mcmc->p.theta->size1, mcmc->p.theta->size2);
  mcmc->w.theta_p      = gsl_matrix_calloc(mcmc->p.theta->size1, mcmc->p.theta->size2);
  mcmc->w.theta_star   = gsl_matrix_calloc(mcmc->p.theta->size1, mcmc->p.theta->size2);
  mcmc->w.theta_star_p = gsl_matrix_calloc(mcmc->p.theta->size1, mcmc->p.theta->size2);
  mcmc->w.theta_tune   = gsl_matrix_calloc(mcmc->p.theta->size1, mcmc->p.theta->size2);

  gsl_matrix_memcpy(mcmc->w.theta,      mcmc->p.theta);
  gsl_matrix_memcpy(mcmc->w.theta_star, mcmc->p.theta);

  gsl_matrix_set_all(mcmc->w.theta_tune, mcmc->p.theta_tune_sd);
}

static
void
qdm_mcmc_workspace_fini(
    qdm_mcmc *mcmc
)
{
  gsl_vector_free(mcmc->w.xi);
  gsl_vector_free(mcmc->w.xi_acc);
  gsl_vector_free(mcmc->w.xi_p);
  gsl_vector_free(mcmc->w.xi_tune);

  gsl_vector_free(mcmc->w.ll);
  gsl_vector_free(mcmc->w.ll_p);

  gsl_vector_free(mcmc->w.tau);
  gsl_vector_free(mcmc->w.tau_p);

  gsl_matrix_free(mcmc->w.theta);
  gsl_matrix_free(mcmc->w.theta_acc);
  gsl_matrix_free(mcmc->w.theta_p);
  gsl_matrix_free(mcmc->w.theta_star);
  gsl_matrix_free(mcmc->w.theta_star_p);
  gsl_matrix_free(mcmc->w.theta_tune);

  mcmc->w.xi      = NULL;
  mcmc->w.xi_acc  = NULL;
  mcmc->w.xi_p    = NULL;
  mcmc->w.xi_tune = NULL;

  mcmc->w.ll   = NULL;
  mcmc->w.ll_p = NULL;

  mcmc->w.tau   = NULL;
  mcmc->w.tau_p = NULL;

  mcmc->w.theta        = NULL;
  mcmc->w.theta_acc    = NULL;
  mcmc->w.theta_p      = NULL;
  mcmc->w.theta_star   = NULL;
  mcmc->w.theta_star_p = NULL;
  mcmc->w.theta_tune   = NULL;
}

static
void
qdm_mcmc_results_init(
    qdm_mcmc *mcmc
)
{
  mcmc->r.s = mcmc->p.iter / mcmc->p.thin;

  mcmc->r.theta      = qdm_ijk_calloc(mcmc->p.theta->size1, mcmc->p.theta->size2, mcmc->r.s);
  mcmc->r.theta_star = qdm_ijk_calloc(mcmc->p.theta->size1, mcmc->p.theta->size2, mcmc->r.s);

  mcmc->r.ll  = qdm_ijk_calloc(1, mcmc->p.ll->size, mcmc->r.s);
  mcmc->r.tau = qdm_ijk_calloc(1, mcmc->p.tau->size, mcmc->r.s);

  mcmc->r.xi = qdm_ijk_calloc(1, mcmc->p.xi->size, mcmc->r.s);
}

static
void
qdm_mcmc_results_fini(
    qdm_mcmc *mcmc
)
{
  qdm_ijk_free(mcmc->r.theta);
  qdm_ijk_free(mcmc->r.theta_star);

  qdm_ijk_free(mcmc->r.ll);
  qdm_ijk_free(mcmc->r.tau);

  qdm_ijk_free(mcmc->r.xi);

  mcmc->r.theta      = NULL;
  mcmc->r.theta_star = NULL;

  mcmc->r.ll  = NULL;
  mcmc->r.tau = NULL;

  mcmc->r.xi = NULL;
}

qdm_mcmc *
qdm_mcmc_alloc(
    qdm_mcmc_parameters p
)
{
  qdm_mcmc *mcmc = malloc(sizeof(qdm_mcmc));

  mcmc->p = p;

  qdm_mcmc_workspace_init(mcmc);
  qdm_mcmc_results_init(mcmc);

  return mcmc;
}

void
qdm_mcmc_free(
    qdm_mcmc *mcmc
)
{
  if (mcmc == NULL) {
    return;
  }

  qdm_mcmc_results_fini(mcmc);
  qdm_mcmc_workspace_fini(mcmc);

  free(mcmc);
}

int
qdm_mcmc_run(
    qdm_mcmc *mcmc
)
{
  int status = 0;

  for (size_t i = 0; i < mcmc->p.burn; i++) {
    status = qdm_mcmc_next(mcmc);
    if (status != 0) {
      goto cleanup;
    }

    if (i % mcmc->p.acc_check == 0) {
      qdm_mcmc_update_tune(mcmc);
    }
  }

  for (size_t i = 0; i < mcmc->p.iter; i++) {
    status = qdm_mcmc_next(mcmc);
    if (status != 0) {
      goto cleanup;
    }

    if (i % mcmc->p.thin == 0) {
      status = qdm_mcmc_save(mcmc, i / mcmc->p.thin);
      if (status != 0) {
        goto cleanup;
      }
    }
  }

cleanup:
  return status;
}

int
qdm_mcmc_next(
    qdm_mcmc *mcmc
)
{
  int status = 0;

  status = qdm_mcmc_update_theta(mcmc);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_mcmc_update_xi(mcmc);
  if (status != 0) {
    goto cleanup;
  }

cleanup:
  return status;
}

void
qdm_mcmc_update_tune(
    qdm_mcmc *mcmc
)
{
  for (size_t p = 0; p < mcmc->w.theta->size1; p++) {
    for (size_t m = 0; m < mcmc->w.theta->size2; m++) {
      if (gsl_matrix_get(mcmc->w.theta_acc, p, m) / (double)mcmc->p.acc_check > 0.5) {
        gsl_matrix_set(mcmc->w.theta_tune, p, m, gsl_min(gsl_matrix_get(mcmc->w.theta_tune, p, m) * 1.2, 5));
      }

      if (gsl_matrix_get(mcmc->w.theta_acc, p, m) / (double)mcmc->p.acc_check < 0.3) {
        gsl_matrix_set(mcmc->w.theta_tune, p, m, gsl_matrix_get(mcmc->w.theta_tune, p, m) * 0.8);
      }
    }
  }

  if (gsl_vector_get(mcmc->w.xi_acc, 0) / mcmc->p.acc_check > 0.5) {
    gsl_vector_set(mcmc->w.xi_tune, 0, gsl_min(gsl_vector_get(mcmc->w.xi_tune, 0) * 1.2, 5));
  }
  if (gsl_vector_get(mcmc->w.xi_acc, 0) / mcmc->p.acc_check < 0.3) {
    gsl_vector_set(mcmc->w.xi_tune, 0, gsl_vector_get(mcmc->w.xi_tune, 0) * 0.8);
  }

  if (gsl_vector_get(mcmc->w.xi_acc, 1) / mcmc->p.acc_check > 0.5) {
    gsl_vector_set(mcmc->w.xi_tune, 1, gsl_min(gsl_vector_get(mcmc->w.xi_tune, 1) * 1.2, 5));
  }
  if (gsl_vector_get(mcmc->w.xi_acc, 1) / mcmc->p.acc_check < 0.3) {
    gsl_vector_set(mcmc->w.xi_tune, 1, gsl_vector_get(mcmc->w.xi_tune, 1) * 0.8);
  }

  /* Reset the acceptance counter. */
  gsl_matrix_set_zero(mcmc->w.theta_acc);
  gsl_vector_set_zero(mcmc->w.xi_acc);
}

int
qdm_mcmc_update_theta(
    qdm_mcmc *mcmc
)
{
  int status = 0;

  /* Update the mth spline. */
  size_t p_max = mcmc->w.theta->size1;
  if (mcmc->p.truncate) {
    p_max = 1;
  }

  for (size_t p = 0; p < p_max; p++) {
    for (size_t m = 0; m < mcmc->w.theta->size2; m++) {
      /* Propose new theta star. */
      status = gsl_matrix_memcpy(mcmc->w.theta_star_p, mcmc->w.theta_star);
      if (status != 0) {
        goto cleanup;
      }

      double proposal =
        gsl_matrix_get(mcmc->w.theta_tune, p, m) *
        gsl_ran_ugaussian(mcmc->p.rng) +
        gsl_matrix_get(mcmc->w.theta_star, p, m);

      gsl_matrix_set(mcmc->w.theta_star_p, p, m, proposal);

      gsl_matrix_memcpy(mcmc->w.theta_p, mcmc->w.theta_star_p);
      qdm_theta_matrix_constrain(mcmc->w.theta_p, mcmc->p.theta_min);

      for (size_t i = 0; i < mcmc->p.y->size; i++) {
        qdm_logl_3(
            gsl_vector_ptr(mcmc->w.ll_p, i),
            gsl_vector_ptr(mcmc->w.tau_p, i),

            gsl_vector_get(mcmc->p.x, i),
            gsl_vector_get(mcmc->p.y, i),

            mcmc->p.t,
            mcmc->w.xi,

            mcmc->w.theta_p
        );
      }

      double ratio = qdm_vector_sum(mcmc->w.ll_p) - qdm_vector_sum(mcmc->w.ll);

      if (log(gsl_rng_uniform(mcmc->p.rng)) < ratio) {
        status = gsl_matrix_memcpy(mcmc->w.theta_star, mcmc->w.theta_star_p);
        if (status != 0) {
          goto cleanup;
        }

        status = gsl_matrix_memcpy(mcmc->w.theta, mcmc->w.theta_p);
        if (status != 0) {
          goto cleanup;
        }

        status = gsl_vector_memcpy(mcmc->w.ll, mcmc->w.ll_p);
        if (status != 0) {
          goto cleanup;
        }

        status = gsl_vector_memcpy(mcmc->w.tau, mcmc->w.tau_p);
        if (status != 0) {
          goto cleanup;
        }

        gsl_matrix_set(mcmc->w.theta_acc, p, m, gsl_matrix_get(mcmc->w.theta_acc, p, m) + 1);
      }
    }
  }

cleanup:
  return status;
}

int
qdm_mcmc_update_xi(
    qdm_mcmc *mcmc
)
{
  int status = 0;

  /* Update xi_low. */
  {
    double xi_low = gsl_vector_get(mcmc->w.xi, 0);
    double xi_low_p = exp(log(xi_low) + gsl_vector_get(mcmc->w.xi_tune, 0) * gsl_ran_ugaussian(mcmc->p.rng));

    gsl_vector_memcpy(mcmc->w.xi_p, mcmc->w.xi);
    gsl_vector_set(mcmc->w.xi_p, 0, xi_low_p);

    for (size_t i = 0; i < mcmc->p.y->size; i++) {
      qdm_logl_3(
          gsl_vector_ptr(mcmc->w.ll_p, i),
          gsl_vector_ptr(mcmc->w.tau_p, i),

          gsl_vector_get(mcmc->p.x, i),
          gsl_vector_get(mcmc->p.y, i),

          mcmc->p.t,
          mcmc->w.xi_p,

          mcmc->w.theta
      );
    }

    double ratio = -0.5 * (1 / mcmc->p.xi_prior_var) *  pow(log(xi_low_p) - mcmc->p.xi_prior_mean, 2) +
                    0.5 * (1 / mcmc->p.xi_prior_var) * (pow(log(xi_low  ) - mcmc->p.xi_prior_mean, 2) + qdm_vector_sum(mcmc->w.ll_p) - qdm_vector_sum(mcmc->w.ll));
    if (log(gsl_rng_uniform(mcmc->p.rng)) < ratio) {
      gsl_vector_set(mcmc->w.xi, 0, xi_low_p);

      status = gsl_vector_memcpy(mcmc->w.ll, mcmc->w.ll_p);
      if (status != 0) {
        goto cleanup;
      }

      gsl_vector_set(mcmc->w.xi_acc, 0, gsl_vector_get(mcmc->w.xi_acc, 0) + 1);
    }
  }

  /* Update xi_high. */
  {
    double xi_high = gsl_vector_get(mcmc->w.xi, 1);
    double xi_high_p = exp(log(xi_high) + gsl_vector_get(mcmc->w.xi_tune, 1) * gsl_ran_ugaussian(mcmc->p.rng));

    gsl_vector_memcpy(mcmc->w.xi_p, mcmc->w.xi);
    gsl_vector_set(mcmc->w.xi_p, 1, xi_high_p);

    for (size_t i = 0; i < mcmc->p.y->size; i++) {
      qdm_logl_3(
          gsl_vector_ptr(mcmc->w.ll_p, i),
          gsl_vector_ptr(mcmc->w.tau_p, i),

          gsl_vector_get(mcmc->p.x, i),
          gsl_vector_get(mcmc->p.y, i),

          mcmc->p.t,
          mcmc->w.xi_p,

          mcmc->w.theta
      );
    }

    double ratio = -0.5 * (1 / mcmc->p.xi_prior_var) *  pow(log(xi_high_p) - mcmc->p.xi_prior_mean, 2) +
                    0.5 * (1 / mcmc->p.xi_prior_var) * (pow(log(xi_high  ) - mcmc->p.xi_prior_mean, 2) + qdm_vector_sum(mcmc->w.ll_p) - qdm_vector_sum(mcmc->w.ll));
    if (log(gsl_rng_uniform(mcmc->p.rng)) < ratio) {
      gsl_vector_set(mcmc->w.xi, 1, xi_high_p);

      status = gsl_vector_memcpy(mcmc->w.ll, mcmc->w.ll_p);
      if (status != 0) {
        goto cleanup;
      }

      gsl_vector_set(mcmc->w.xi_acc, 1, gsl_vector_get(mcmc->w.xi_acc, 1) + 1);
    }
  }

cleanup:
  return status;
}

int
qdm_mcmc_save(
    qdm_mcmc *mcmc,
    size_t k
)
{
  int status = 0;

  {
    gsl_matrix_view m = qdm_ijk_get_ij(mcmc->r.theta, k);

    status = gsl_matrix_memcpy(&m.matrix, mcmc->w.theta);
    if (status != 0) {
      goto cleanup;
    }
  }

  {
    gsl_matrix_view m = qdm_ijk_get_ij(mcmc->r.theta_star, k);

    status = gsl_matrix_memcpy(&m.matrix, mcmc->w.theta_star);
    if (status != 0) {
      goto cleanup;
    }
  }

  {
    gsl_matrix_view m = qdm_ijk_get_ij(mcmc->r.ll, k);

    status = gsl_matrix_set_row(&m.matrix, 0, mcmc->w.ll);
    if (status != 0) {
      goto cleanup;
    }
  }

  {
    gsl_matrix_view m = qdm_ijk_get_ij(mcmc->r.tau, k);

    status = gsl_matrix_set_row(&m.matrix, 0, mcmc->w.tau);
    if (status != 0) {
      goto cleanup;
    }
  }

  {
    gsl_matrix_view m = qdm_ijk_get_ij(mcmc->r.xi, k);

    status = gsl_matrix_set_row(&m.matrix, 0, mcmc->w.xi);
    if (status != 0) {
      goto cleanup;
    }
  }

cleanup:
  return status;
}

int
qdm_mcmc_write(
    hid_t id,
    const qdm_mcmc *mcmc
)
{
  int status = 0;

  status = qdm_ijk_write(id, "theta", mcmc->r.theta);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_ijk_write(id, "theta_star", mcmc->r.theta_star);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_ijk_write(id, "ll", mcmc->r.ll);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_ijk_write(id, "tau", mcmc->r.tau);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_ijk_write(id, "xi", mcmc->r.xi);
  if (status != 0) {
    goto cleanup;
  }

cleanup:
  return status;
}

int
qdm_mcmc_read(
    hid_t id,
    qdm_mcmc **mcmc
)
{
  int status = 0;

  *mcmc = malloc(sizeof(qdm_mcmc));

  status = qdm_ijk_read(id, "theta", &(*mcmc)->r.theta);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_ijk_read(id, "theta_star", &(*mcmc)->r.theta_star);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_ijk_read(id, "ll", &(*mcmc)->r.ll);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_ijk_read(id, "tau", &(*mcmc)->r.tau);
  if (status != 0) {
    goto cleanup;
  }

  status = qdm_ijk_read(id, "xi", &(*mcmc)->r.xi);
  if (status != 0) {
    goto cleanup;
  }

  (*mcmc)->r.s = (*mcmc)->r.theta->size3;

cleanup:
  return status;
}
