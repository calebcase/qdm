#include "munit.h"

#include "test_gsl.c"
#include "test_knots.c"

static MunitSuite test_suites[] = {
  test_suite_gsl,
  test_suite_knots,
  { NULL, NULL, NULL, 0, MUNIT_SUITE_OPTION_NONE },
};

static const MunitSuite test_suite = {
  ""                      , // name
  NULL                    , // tests
  test_suites             , // suites
  1                       , // iterations
  MUNIT_SUITE_OPTION_NONE , // options
};

int
main(int argc, char** argv)
{
  return munit_suite_main(&test_suite, NULL, argc, argv);
}
