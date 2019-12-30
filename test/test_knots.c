#include <gsl/gsl_sort_vector.h>
#include <munit.h>

#include <qdm.h>

#include "test.h"

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

  for (size_t i = 0; i < 3; i++) {
    munit_assert_double(gsl_vector_get(knots, i), ==, 0);
  }

  for (size_t i = spline_df + knots_inter->size; i < spline_df * 2 + knots_inter->size; i++) {
    munit_assert_double(gsl_vector_get(knots, i), ==, 1);
  }

  return MUNIT_OK;
}

static
MunitResult
test_qdm_knots_optimize(
    const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  int status = 0;

  double tau_low = 0.01;
  double tau_high = 0.99;

  gsl_vector *result = gsl_vector_calloc(2);

  const char *y_path = munit_parameters_get(params, "y_path");
  munit_logf(MUNIT_LOG_DEBUG, "Y Path: %s", y_path);

  FILE *y_file = fopen(y_path, "r");
  munit_assert_ptr(y_file, !=, NULL);

  gsl_vector *sorted_data = qdm_vector_csv_fread(y_file);
  gsl_sort_vector(sorted_data);
  fclose(y_file);

  gsl_vector *middle = qdm_vector_seq(tau_low, tau_high, 0.005);

  double possible_low = tau_low + 0.1;
  double possible_high = tau_high - 0.1;
  double possible_delta = (possible_high - possible_low) / 19;

  gsl_vector *possible_knots = qdm_vector_seq(possible_low, possible_high, possible_delta);
  munit_log_vector(MUNIT_LOG_DEBUG, possible_knots);

  size_t iterate_n = 1000;
  size_t spline_df = 3;

  munit_log_vector_dim(MUNIT_LOG_DEBUG, result);
  munit_log_vector_dim(MUNIT_LOG_DEBUG, sorted_data);
  munit_log_vector_dim(MUNIT_LOG_DEBUG, middle);
  munit_log_vector_dim(MUNIT_LOG_DEBUG, possible_knots);
  munit_log_size(MUNIT_LOG_DEBUG, iterate_n);
  munit_log_size(MUNIT_LOG_DEBUG, spline_df);

  status = qdm_knots_optimize(
      result,
      sorted_data,
      middle,
      possible_knots,

      iterate_n,
      spline_df
  );
  munit_assert_int(status, ==, 0);

  munit_log_vector(MUNIT_LOG_DEBUG, result);

  munit_assert_double_equal(gsl_vector_get(result, 0), 0.15105263157894738, 10);
  munit_assert_double_equal(gsl_vector_get(result, 1), 0.89000000000000001, 10);

  return MUNIT_OK;
}

static MunitParameterEnum test_params_qdm_knots_optimize[] = {
  { "y_path" , NULL } ,
  { NULL     , NULL } ,
};

MunitTest tests_knots[] = {
  {
    "/qdm_knots_vector"    , // name
    test_qdm_knots_vector  , // test
    NULL                   , // setup
    NULL                   , // tear_down
    MUNIT_TEST_OPTION_NONE , // options
    NULL                   , // parameters
  },
  {
    "/qdm_knots_optimize"          , // name
    test_qdm_knots_optimize        , // test
    NULL                           , // setup
    NULL                           , // tear_down
    MUNIT_TEST_OPTION_NONE         , // options
    test_params_qdm_knots_optimize , // parameters
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
