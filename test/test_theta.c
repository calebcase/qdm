#include <munit.h>

#include <qdm.h>

#include "test.h"

static
MunitResult
test_qdm_theta_matrix_constrain(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  double min = 1e-04;

  /* Basic (unmodified) */
  double basic_data[] = {
    11.2382694454771581,  4.89597270237975835,  7.53285823736484961,  0.0449150095737098978, 3.75896703535460786, 1.86638774220197923,
     1.5098447359663183, -0.42943079543658147, -0.13782001356078322, -0.0091277541152237976, 0.55811649692923071, 0.28116069158352697,
  };
  gsl_matrix_view basic = gsl_matrix_view_array(basic_data, 2, 6);

  qdm_theta_matrix_constrain(&basic.matrix, min);

  double expected_basic_data[] = {
    11.2382694454771581,  4.89597270237975835,  7.53285823736484961,  0.0449150095737098978, 3.75896703535460786, 1.86638774220197923,
     1.5098447359663183, -0.42943079543658147, -0.13782001356078322, -0.0091277541152237976, 0.55811649692923071, 0.28116069158352697,
  };
  gsl_matrix_view expected_basic = gsl_matrix_view_array(expected_basic_data, 2, 6);

  munit_matrix_equal(&basic.matrix, &expected_basic.matrix, 10);

  /* Zero */
  double zero_data[] = {
      1,   2,   3,   4,   5,   6,
    -10, -20, -30, -40, -50, -60,
  };
  gsl_matrix_view zero = gsl_matrix_view_array(zero_data, 2, 6);

  qdm_theta_matrix_constrain(&zero.matrix, min);

  double expected_zero_data[] = {
      1, 2, 3, 4, 5, 6,
    -10, 0, 0, 0, 0, 0,
  };
  gsl_matrix_view expected_zero = gsl_matrix_view_array(expected_zero_data, 2, 6);

  munit_matrix_equal(&zero.matrix, &expected_zero.matrix, 10);

  /* Zero More */
  double zero_more_data[] = {
      1,  2,  3,  4,  5,  6,
     -1, -1, -1, -1, -1, -1,
     -1, -1, -1, -1, -1, -1,
     -1, -1, -1, -1, -1, -1,
     -1, -1, -1, -1, -1, -1,
     -1, -1, -1, -1, -1, -1,
     -1, -1, -1, -1, -1, -1,
     -1, -1, -1, -1, -1, -1,
  };
  gsl_matrix_view zero_more = gsl_matrix_view_array(zero_more_data, 8, 6);

  qdm_theta_matrix_constrain(&zero_more.matrix, min);

  double expected_zero_more_data[] = {
      1,  2,  3,  4,  5,  6,
     -1,  0,  0,  0,  0,  0,
     -1,  0,  0,  0,  0, -1,
     -1,  0,  0,  0, -1, -1,
     -1,  0,  0, -1, -1, -1,
     -1,  0, -1, -1, -1, -1,
     -1, -1, -1, -1, -1, -1,
     -1, -1, -1, -1, -1, -1,
  };
  gsl_matrix_view expected_zero_more = gsl_matrix_view_array(expected_zero_more_data, 8, 6);

  munit_matrix_equal(&zero_more.matrix, &expected_zero_more.matrix, 10);

  return MUNIT_OK;
}

MunitTest tests_theta[] = {
  {
    "/qdm_theta_matrix_constrain"   , // name
    test_qdm_theta_matrix_constrain , // test
    NULL                            , // setup
    NULL                            , // tear_down
    MUNIT_TEST_OPTION_NONE          , // options
    NULL                            , // parameters
  },
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL },
};

static const MunitSuite test_suite_theta = {
  "/theta"                , // name
  tests_theta             , // tests
  NULL                    , // suites
  1                       , // iterations
  MUNIT_SUITE_OPTION_NONE , // options
};
