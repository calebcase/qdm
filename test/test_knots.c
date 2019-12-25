#include <munit.h>
#include <qdm.h>

static
MunitResult
test_qdm_knots_vector(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  size_t spline_df = 3;
  gsl_vector *knots_inter = qdm_vector_seq(0.1, 0.9, 0.2);

  gsl_vector *knots = qdm_knots_vector(spline_df, knots_inter);

  munit_assert_size(spline_df * 2 + knots_inter->size, ==, knots->size);

  return MUNIT_OK;
}

MunitTest tests_knots[] = {
  {
    "/qdm_knots_vector"    , // name
    test_qdm_knots_vector  , // test
    NULL                   , // setup
    NULL                   , // tear_down
    MUNIT_TEST_OPTION_NONE , // options
    NULL                   , // parameters
  },
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL },
};

static const MunitSuite test_suite_knots = {
  "/knots"                , // name
  tests_knots             , // tests
  NULL                    , // suites
  1                       , // iterations
  MUNIT_SUITE_OPTION_NONE , // options
};
