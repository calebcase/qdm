#include "munit.h"

#include "test_gsl.c"
#include "test_ijk.c"
#include "test_ispline.c"
#include "test_knots.c"
#include "test_logl.c"
#include "test_mspline.c"
#include "test_tau.c"
#include "test_theta.c"

static MunitSuite test_suites[] = {
  test_suite_gsl,
  test_suite_ijk,
  test_suite_ispline,
  test_suite_knots,
  test_suite_logl,
  test_suite_mspline,
  test_suite_tau,
  test_suite_theta,
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
