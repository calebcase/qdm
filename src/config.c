#include <inttypes.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <inih/ini.h>

typedef struct {
#define CFG(s, n, t, d) t s##_##n;
#include "config.def"
#undef CFG
} qdm_config;

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

static void _dump_config_str    ( char *x  )  { printf ( "%s"    , x ) ; }
static void _dump_config_size   ( size_t x )  { printf ( "%zu"   , x ) ; }
static void _dump_config_double ( double x )  { printf ( "%.17g" , x ) ; }

static
void
dump_config(qdm_config *cfg)
{
#define CFG(s, n, t, d) { \
  printf("%s_%s = ", #s, #n); \
  _Generic(((t)d), \
    char *: _dump_config_str, \
    size_t: _dump_config_size, \
    double: _dump_config_double \
  )(cfg->s##_##n); \
  printf("\n"); \
}
#include "config.def"
#undef CFG
}

int
qdm_config_init()
{
  int status = 0;

  qdm_config cfg = {
#define CFG(s, n, t, d) d,
#include "config.def"
#undef CFG
  };

  status = ini_parse("test.ini", qdm_config_handler, &cfg);
  dump_config(&cfg);

  return status;
}
