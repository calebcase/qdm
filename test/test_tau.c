#include <munit.h>

#include <qdm.h>

#include "test.h"

static
MunitResult
test_qdm_tau_find(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  double tau_low = 0.01;
  double tau_high = 0.98999999999999999;
  size_t spline_df = 3;

  double knots_data[] = {
      0, 0, 0,
      0.15105263157894738,
      0.89000000000000001,
      1, 1, 1,
  };
  gsl_vector_view knots = gsl_vector_view_array(knots_data, sizeof(knots_data) / sizeof(double));
  munit_log_vector_dim(MUNIT_LOG_DEBUG, &knots.vector);

  qdm_tau *t = qdm_tau_alloc(tau_low, tau_high, spline_df, &knots.vector);

  double reset_data[] = {
    0.010000000000000000,
    0.010100000000000000,
    0.010200000000000001,
    0.010300000000000000,
    0.010400000000000000,
    0.010500000000000001,
    0.010600000000000000,
    0.010699999999999999,
    0.010800000000000001,
    0.010900000000000000,
    0.010999999999999999,
    0.010999999999999999,
    0.150714285714285717,
    0.290428571428571425,
    0.430142857142857160,
    0.569857142857142840,
    0.709571428571428520,
    0.849285714285714310,
    0.988999999999999990,
    0.988999999999999990,
    0.989099999999999979,
    0.989199999999999968,
    0.989299999999999957,
    0.989399999999999946,
    0.989500000000000046,
    0.989600000000000035,
    0.989700000000000024,
    0.989800000000000013,
    0.989900000000000002,
    0.989999999999999991,
  };
  gsl_vector_view reset = gsl_vector_view_array(reset_data, sizeof(reset_data) / sizeof(double));
  munit_vector_equal(&reset.vector, t->reset, 3);

  double v = 13.2690000000000001;

  double mmm_data[] = {
      11.238269445477158115,
       4.895972702379758346,
       7.532858237364849607,
       0.044915009573709898,
       3.758967035354607855,
       1.866387742201979227,
  };
  gsl_vector_view mmm = gsl_vector_view_array(mmm_data, sizeof(mmm_data) / sizeof(double));
  munit_log_vector_dim(MUNIT_LOG_DEBUG, &mmm.vector);

  double result = qdm_tau_find(t, v, &mmm.vector);

  double expected = 0.023441441362023058;

  munit_assert_double_equal(expected, result, 5);

  qdm_tau_free(t);

  return MUNIT_OK;
}

MunitTest tests_tau[] = {
  {
    "/qdm_tau_find",        // name
    test_qdm_tau_find,      // test
    NULL,                   // setup
    NULL,                   // tear_down
    MUNIT_TEST_OPTION_NONE, // options
    NULL,                   // parameters
  },
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL },
};

static const MunitSuite test_suite_tau = {
  "/tau",                  // name
  tests_tau,               // tests
  NULL,                    // suites
  1,                       // iterations
  MUNIT_SUITE_OPTION_NONE, // options
};
