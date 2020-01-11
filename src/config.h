#ifndef QDM_CONFIG_H
#define QDM_CONFIG_H 1

typedef struct {
#define CFG(s, n, t, d) t s##_##n;
#include "config.def"
#undef CFG
} qdm_config;

void
qdm_config_fprint(FILE *f, qdm_config *cfg);

int
qdm_config_new(FILE *f, qdm_config *cfg);

#endif /* QDM_CONFIG_H */
