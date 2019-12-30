#include <munit.h>

#include <qdm.h>

#include "test.h"

static
MunitResult
test_qdm_mspline_vector(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  double tau = 0.014999999999999999;
  size_t spline_df = 3;

  double knots_data[] = {
      0, 0, 0,
      0.31526315789473686,
      0.47947368421052627,
      1, 1, 1,
  };
  gsl_vector_view knots = gsl_vector_view_array(knots_data, sizeof(knots_data) / sizeof(double));
  munit_log_vector_dim(MUNIT_LOG_DEBUG, &knots.vector);

  gsl_vector *result = gsl_vector_alloc((knots.vector.size - spline_df) + 1);
  munit_log_vector_dim(MUNIT_LOG_DEBUG, result);

  double expected_data[] = {
      0,
      8.6318858004720109,
      0.57191661896863488,
      0.0044654555983353154,
      0,
      0,
  };
  gsl_vector_view expected = gsl_vector_view_array(expected_data, sizeof(expected_data) / sizeof(double));
  munit_log_vector_dim(MUNIT_LOG_DEBUG, &expected.vector);

  munit_assert_size(expected.vector.size, ==, result->size);

  qdm_mspline_vector(result, tau, spline_df, &knots.vector);
  munit_log_vector(MUNIT_LOG_DEBUG, result);

  for (size_t i = 0; i < expected.vector.size; i++) {
    double e = gsl_vector_get(&expected.vector, i);
    double r = gsl_vector_get(result, i);

    munit_logf(MUNIT_LOG_DEBUG, "(%zu) %0.17g ?= %0.17g (delta = %0.17g)", i, e, r, fabs(e - r));
    munit_assert_double_equal(e, r, 10);
  }

  return MUNIT_OK;
}

MunitTest tests_mspline[] = {
  {
    "/qdm_mspline_vector"   , // name
    test_qdm_mspline_vector , // test
    NULL                    , // setup
    NULL                    , // tear_down
    MUNIT_TEST_OPTION_NONE  , // options
    NULL                    , // parameters
  },
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL },
};

static const MunitSuite test_suite_mspline = {
  "/mspline"              , // name
  tests_mspline           , // tests
  NULL                    , // suites
  1                       , // iterations
  MUNIT_SUITE_OPTION_NONE , // options
};
