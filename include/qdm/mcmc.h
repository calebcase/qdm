#ifndef QDM_MCMC_H
#define QDM_MCMC_H 1

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <hdf5.h>

#include "ijk.h"
#include "tau.h"

typedef struct {
  const gsl_rng *rng;

  size_t burn;
  size_t iter;
  size_t thin;

  size_t acc_check;
  size_t spline_df;

  const gsl_vector *xi;
  const gsl_vector *ll;
  const gsl_vector *tau;
  const gsl_matrix *theta;

  double xi_prior_mean;
  double xi_prior_var;
  double xi_tune_sd;

  double theta_min;
  double theta_tune_sd;

  const gsl_vector *x;
  const gsl_vector *y;

  const qdm_tau *t;
  const gsl_vector *m_knots;

  bool truncate;
} qdm_mcmc_parameters;

typedef struct {
  gsl_vector *xi;
  gsl_vector *xi_acc;
  gsl_vector *xi_p;
  gsl_vector *xi_tune;

  gsl_vector *ll;
  gsl_vector *ll_p;

  gsl_vector *tau;
  gsl_vector *tau_p;

  gsl_matrix *theta;
  gsl_matrix *theta_acc;
  gsl_matrix *theta_p;
  gsl_matrix *theta_star;
  gsl_matrix *theta_star_p;
  gsl_matrix *theta_tune;
} qdm_mcmc_workspace;

typedef struct {
  size_t s;

  qdm_ijk *theta;
  qdm_ijk *theta_star;

  qdm_ijk *ll;
  qdm_ijk *tau;

  qdm_ijk *xi;
} qdm_mcmc_results;

typedef struct {
  qdm_mcmc_parameters p;
  qdm_mcmc_workspace w;
  qdm_mcmc_results r;
} qdm_mcmc;

qdm_mcmc *
qdm_mcmc_alloc(
    qdm_mcmc_parameters p
);

void
qdm_mcmc_free(
    qdm_mcmc *mcmc
);

int
qdm_mcmc_run(
    qdm_mcmc *mcmc
);

int
qdm_mcmc_next(
    qdm_mcmc *mcmc
);

int
qdm_mcmc_update_theta(
    qdm_mcmc *mcmc
);

int
qdm_mcmc_update_xi(
    qdm_mcmc *mcmc
);

void
qdm_mcmc_update_tune(
    qdm_mcmc *mcmc
);

int
qdm_mcmc_save(
    qdm_mcmc *mcmc,
    size_t k
);

int
qdm_mcmc_write(
    hid_t id,
    const qdm_mcmc *mcmc
);

int
qdm_mcmc_read(
    hid_t id,
    qdm_mcmc **mcmc
);

#endif /* QDM_MCMC_H */
