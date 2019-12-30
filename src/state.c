#include <gsl/gsl_rng.h>
#include <qdm.h>

qdm_state_t *qdm_state;

static void qdm_state_init() __attribute__((constructor));

static
void
qdm_state_init()
{
  qdm_state = malloc(sizeof(qdm_state_t));

  gsl_rng_env_setup();
  qdm_state->rng = gsl_rng_alloc(gsl_rng_default);
}
