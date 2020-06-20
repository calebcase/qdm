#ifndef QDM_PHASE_H
#define QDM_PHASE_H 1

typedef struct {
  /*     .---.
   *    /   /|
   *   +---+ .
   * 2 |   |/ s
   *   +---+
   *     m
   */
  qdm_ijk *theta;
  qdm_ijk *theta_star;

  /*   +---+
   * 2 |   |
   *   +---+
   *     m
   */
  gsl_matrix *theta_acc;

  /*
   *     .-----------.
   *    /     /     /|
   *   +-----------+ .
   * 1 | low | hi  |/ s
   *   +-----------+
   *         2
   */
  qdm_ijk *xi;

  /*   +-----------+
   * 1 | low | hi  |
   *   +-----------+
   *         2
   */
  gsl_matrix *xi_acc;
} qdm_phase;

qdm_phase *
qdm_phase_alloc();

void
qdm_phase_free(
    qdm_phase *p
);

int
qdm_phase_run(
    qdm_phase *p
);

int
qdm_phase_write(
    const qdm_phase *p,
    hid_t id,
    const char *name
);

#endif /* QDM_PHASE_H */
