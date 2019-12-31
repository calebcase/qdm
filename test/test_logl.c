#include <munit.h>

#include <qdm.h>

#include "test.h"

static
MunitResult
test_qdm_find_tau(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  int status = 0;

  double v = 13.2690000000000001;
  size_t spline_df = 3;

  double knots_data[] = {
      0, 0, 0,
      0.15105263157894738,
      0.89000000000000001,
      1, 1, 1,
  };
  gsl_vector_view knots = gsl_vector_view_array(knots_data, sizeof(knots_data) / sizeof(double));
  munit_log_vector_dim(MUNIT_LOG_DEBUG, &knots.vector);

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

  double result = 0;
  double expected = 0.023441441362023058;

  status = qdm_find_tau(&result, v, spline_df, &knots.vector, &mmm.vector);
  munit_assert_int(status, ==, 0);

  munit_assert_double_equal(expected, result, 10);

  return MUNIT_OK;
}

static
MunitResult
test_qdm_logl_less_than(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  int status = 0;

  size_t spline_df = 3;

  double knots_data[] = {
      0, 0, 0,
      0.15105263157894738,
      0.89000000000000001,
      1, 1, 1,
  };
  gsl_vector_view knots = gsl_vector_view_array(knots_data, sizeof(knots_data) / sizeof(double));
  munit_log_vector_dim(MUNIT_LOG_DEBUG, &knots.vector);

  double min = 1e-04;

  double theta_data[] = {
      11.2382694454771581,  4.89597270237975835,  7.53285823736484961,  0.0449150095737098978, 3.75896703535460786, 1.86638774220197923,
       1.5098447359663183, -0.42943079543658147, -0.13782001356078322, -0.0091277541152237976, 0.55811649692923071, 0.28116069158352697,
  };
  gsl_matrix_view theta = gsl_matrix_view_array(theta_data, 2, 6);

  qdm_theta_matrix_constrain(&theta.matrix, min);

  double mmm_data[] = {
      11.238269445477158115,
       4.895972702379758346,
       7.532858237364849607,
       0.044915009573709898,
       3.758967035354607855,
       1.866387742201979227,
  };
  gsl_vector_view mmm = gsl_vector_view_array(mmm_data, 6);

  double v = 10.983000000000001;

  double tau_low = 0.01;
  double tau_high = 0.98999999999999999;

  double xi_low = 0.20000000000000001;
  double xi_high = 0.20000000000000001;

  double log_likelihood = 0;
  double tau = 0;

  status = qdm_logl(
      &log_likelihood,
      &tau,

      v,
      spline_df,
      &knots.vector,
      &mmm.vector,

      tau_low,
      tau_high,
      xi_low,
      xi_high
  );
  munit_assert_int(status, ==, 0);

  double expected_log_likelihood = -5.9039234039068802;
  double expected_tau = 0.0030463905316060647;

  munit_assert_double_equal(log_likelihood, expected_log_likelihood, 10);
  munit_assert_double_equal(tau, expected_tau, 10);

  return MUNIT_OK;
}

static
MunitResult
test_qdm_logl_greater_than(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  int status = 0;

  size_t spline_df = 3;

  double knots_data[] = {
      0, 0, 0,
      0.15105263157894738,
      0.89000000000000001,
      1, 1, 1,
  };
  gsl_vector_view knots = gsl_vector_view_array(knots_data, sizeof(knots_data) / sizeof(double));
  munit_log_vector_dim(MUNIT_LOG_DEBUG, &knots.vector);

  double min = 1e-04;

  double theta_data[] = {
      11.2382694454771581,  4.89597270237975835,  7.53285823736484961,  0.0449150095737098978, 3.75896703535460786, 1.86638774220197923,
       1.5098447359663183, -0.42943079543658147, -0.13782001356078322, -0.0091277541152237976, 0.55811649692923071, 0.28116069158352697,
  };
  gsl_matrix_view theta = gsl_matrix_view_array(theta_data, 2, 6);

  qdm_theta_matrix_constrain(&theta.matrix, min);

  double mmm_data[] = {
      11.823719445137566453,
       4.729458720475777866,
       7.479417823943321331,
       0.041375676345357812,
       3.975379554572064489,
       1.975409234856816187,
  };
  gsl_vector_view mmm = gsl_vector_view_array(mmm_data, 6);

  double v = 30.0829999999999984;

  double tau_low = 0.01;
  double tau_high = 0.98999999999999999;

  double xi_low = 0.20000000000000001;
  double xi_high = 0.20000000000000001;

  double log_likelihood = 0;
  double tau = 0;

  status = qdm_logl(
      &log_likelihood,
      &tau,

      v,
      spline_df,
      &knots.vector,
      &mmm.vector,

      tau_low,
      tau_high,
      xi_low,
      xi_high
  );
  munit_assert_int(status, ==, 0);

  double expected_log_likelihood = -5.1365196361984982;
  double expected_tau = 0.9965798555939035;

  munit_assert_double_equal(log_likelihood, expected_log_likelihood, 10);
  munit_assert_double_equal(tau, expected_tau, 10);

  return MUNIT_OK;
}

static
MunitResult
test_qdm_logl_within(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  int status = 0;

  size_t spline_df = 3;

  double knots_data[] = {
      0, 0, 0,
      0.15105263157894738,
      0.89000000000000001,
      1, 1, 1,
  };
  gsl_vector_view knots = gsl_vector_view_array(knots_data, sizeof(knots_data) / sizeof(double));
  munit_log_vector_dim(MUNIT_LOG_DEBUG, &knots.vector);

  double min = 1e-04;

  double theta_data[] = {
      11.2382694454771581,  4.89597270237975835,  7.53285823736484961,  0.0449150095737098978, 3.75896703535460786, 1.86638774220197923,
       1.5098447359663183, -0.42943079543658147, -0.13782001356078322, -0.0091277541152237976, 0.55811649692923071, 0.28116069158352697,
  };
  gsl_matrix_view theta = gsl_matrix_view_array(theta_data, 2, 6);

  qdm_theta_matrix_constrain(&theta.matrix, min);

  double mmm_data[] = {
      11.238269445477158115,
       4.895972702379758346,
       7.532858237364849607,
       0.044915009573709898,
       3.758967035354607855,
       1.866387742201979227,
  };
  gsl_vector_view mmm = gsl_vector_view_array(mmm_data, 6);

  double v = 13.2690000000000001;

  double tau_low = 0.01;
  double tau_high = 0.98999999999999999;

  double xi_low = 0.20000000000000001;
  double xi_high = 0.20000000000000001;

  double log_likelihood = 0;
  double tau = 0;

  status = qdm_logl(
      &log_likelihood,
      &tau,

      v,
      spline_df,
      &knots.vector,
      &mmm.vector,

      tau_low,
      tau_high,
      xi_low,
      xi_high
  );
  munit_assert_int(status, ==, 0);

  double expected_log_likelihood = -4.3381416615763744;
  double expected_tau = 0.023441441362023058;

  munit_assert_double_equal(log_likelihood, expected_log_likelihood, 10);
  munit_assert_double_equal(tau, expected_tau, 10);

  return MUNIT_OK;
}


MunitTest tests_logl[] = {
  {
    "/qdm_find_tau"        , // name
    test_qdm_find_tau      , // test
    NULL                   , // setup
    NULL                   , // tear_down
    MUNIT_TEST_OPTION_NONE , // options
    NULL                   , // parameters
  },
  {
    "/qdm_logl/less_than"   , // name
    test_qdm_logl_less_than , // test
    NULL                    , // setup
    NULL                    , // tear_down
    MUNIT_TEST_OPTION_NONE  , // options
    NULL                    , // parameters
  },
  {
    "/qdm_logl/greater_than"   , // name
    test_qdm_logl_greater_than , // test
    NULL                       , // setup
    NULL                       , // tear_down
    MUNIT_TEST_OPTION_NONE     , // options
    NULL                       , // parameters
  },
  {
    "/qdm_logl/within"     , // name
    test_qdm_logl_within   , // test
    NULL                   , // setup
    NULL                   , // tear_down
    MUNIT_TEST_OPTION_NONE , // options
    NULL                   , // parameters
  },
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL },
};

static const MunitSuite test_suite_logl = {
  "/logl"                 , // name
  tests_logl              , // tests
  NULL                    , // suites
  1                       , // iterations
  MUNIT_SUITE_OPTION_NONE , // options
};
