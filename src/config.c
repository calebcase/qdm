#include <inttypes.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <inih/ini.h>

#include "config.h"

int
qdm_config_handler(
    void *user,
    const char *section,
    const char *name,
    const char *value
)
{
  qdm_config *cfg = (qdm_config *)user;

#define CMP(s, n) if (strcmp(section, #s) == 0 && strcmp(name, #n) == 0)
#define CFG(s, n, t, d) CMP(s, n) { \
  cfg->s##_##n = \
    _Generic(((t)d), \
      char *: strdup(value), \
      size_t: strtoumax(value, NULL, 10), \
      double: strtod(value, NULL) \
    ); \
  return 1; \
}
#include "config.def"
#undef CFG
#undef CMP

  return 0;
}

static void _qdm_config_fprint_str   (FILE *f, char *x ) { fprintf(f, "%s"   , x); }
static void _qdm_config_fprint_size  (FILE *f, size_t x) { fprintf(f, "%zu"  , x); }
static void _qdm_config_fprint_double(FILE *f, double x) { fprintf(f, "%.17g", x); }

void
qdm_config_fprint(FILE *f, qdm_config *cfg)
{
#define CFG(s, n, t, d) { \
  fprintf(f, "%s_%s = ", #s, #n); \
  _Generic(((t)d), \
    char *: _qdm_config_fprint_str, \
    size_t: _qdm_config_fprint_size, \
    double: _qdm_config_fprint_double \
  )(f, cfg->s##_##n); \
  fprintf(f, "\n"); \
}
#include "config.def"
#undef CFG
}

int
qdm_config_new(FILE *f, qdm_config *cfg)
{
  qdm_config defaults = {
#define CFG(s, n, t, d) d,
#include "config.def"
#undef CFG
  };
  *cfg = defaults;

  return ini_parse_file(f, qdm_config_handler, cfg);
}
