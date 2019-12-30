#ifndef QDM_STATE_H
#define QDM_STATE_H 1

#include <gsl/gsl_rng.h>

typedef struct {
  gsl_rng *rng;
} qdm_state_t;

extern qdm_state_t *qdm_state;

#endif /* QDM_STATE_H */
