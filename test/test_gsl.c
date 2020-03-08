#include <gsl/gsl_spmatrix.h>
#include <munit.h>
#include <osqp/cs.h>
#include <osqp/osqp.h>
#include <unistd.h>

#include <qdm.h>

#include "test.h"

static
MunitResult
test_qdm_vector_seq(
    const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  double min = atof(munit_parameters_get(params, "qdm_vector_seq_min"));
  double max = atof(munit_parameters_get(params, "qdm_vector_seq_max"));
  double delta = atof(munit_parameters_get(params, "qdm_vector_seq_delta"));

  if (min >= max) {
    return MUNIT_SKIP;
  }

  gsl_vector *seq = qdm_vector_seq(min, max, delta);

  munit_assert_size(((max - min) / delta) + 1, ==, seq->size);

  for (size_t i = 0; i < seq->size; i++) {
    double expecting = min + delta * i;
    double found = gsl_vector_get(seq, i);

    munit_logf(MUNIT_LOG_DEBUG, "%.17g", found);
    munit_assert_double_equal(expecting, found, 3);
  }

  return MUNIT_OK;
}

static
MunitResult
test_qdm_matrix_det_tmm(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  int status = 0;

  /* Zero */
  double data_zero[4][3] = {
      {1,1,1},
      {2,2,2},
      {3,3,3},
      {4,4,4},
  };

  gsl_matrix_view m = gsl_matrix_view_array((double *)data_zero, 4, 3);

  double det_zero = 0;
  status = qdm_matrix_det_tmm(&m.matrix, &det_zero);
  munit_assert_int(status, ==, 0);
  munit_assert_double_equal(0, det_zero, 3);

  /* Non-zero */
  double data_nonzero[4][3] = {
      {1,2,3},
      {4,5,2},
      {1,1,3},
      {1,1,4},
  };

  m = gsl_matrix_view_array((double *)data_nonzero, 4, 3);

  double det_nonzero = 0;
  status = qdm_matrix_det_tmm(&m.matrix, &det_nonzero);
  munit_assert_int(status, ==, 0);
  munit_assert_double_equal(271, det_nonzero, 3);

  return MUNIT_OK;
}

static
MunitResult
test_qdm_vector_quantile(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  /* Basic */
  double data[10] = {
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
  };

  gsl_vector_view v = gsl_vector_view_array((double *)data, 10);
  gsl_vector *probs = qdm_vector_seq(0, 1, 0.25);
  gsl_vector *quantiles = qdm_vector_quantile(&v.vector, probs);

  munit_log(MUNIT_LOG_DEBUG, "basic");
  munit_log_vector(MUNIT_LOG_DEBUG, &v.vector);
  munit_log_vector(MUNIT_LOG_DEBUG, probs);
  munit_log_vector(MUNIT_LOG_DEBUG, quantiles);
  munit_assert_int(quantiles->size, ==, 5);

  double expected[5] = {
    1,
    3.25,
    5.5,
    7.75,
    10,
  };

  for (size_t i = 0; i < 5; i++) {
    munit_assert_double_equal(expected[i], gsl_vector_get(quantiles, i), 3);
  }

  /* Random */
  double rnorm_data[10] = {
    -1.7792761,
    -0.9436401,
    -0.8845596,
    -0.4213213,
    -0.1831428,
    -0.1497455,
    0.1402677,
    0.4509039,
    0.7415657,
    0.8864738,
  };

  v = gsl_vector_view_array((double *)rnorm_data, 10);
  quantiles = qdm_vector_quantile(&v.vector, probs);

  munit_log(MUNIT_LOG_DEBUG, "random");
  munit_log_vector(MUNIT_LOG_DEBUG, &v.vector);
  munit_log_vector(MUNIT_LOG_DEBUG, probs);
  munit_log_vector(MUNIT_LOG_DEBUG, quantiles);
  munit_assert_int(quantiles->size, ==, 5);

  double rnorm_expected[5] = {
    -1.7792761,
    -0.7687500,
    -0.1664441,
    0.3732449,
    0.8864738,
  };

  for (size_t i = 0; i < 5; i++) {
    munit_assert_double_equal(rnorm_expected[i], gsl_vector_get(quantiles, i), 3);
  }

  return MUNIT_OK;
}

static
MunitResult
test_qdm_vector_sum(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  gsl_vector_view v;
  double sum = 0.0;

  /* Basic */
  double basic_data[3] = {
    1,
    2,
    3,
  };

  v = gsl_vector_view_array((double *)basic_data, 3);
  sum = qdm_vector_sum(&v.vector);
  munit_assert_double(6, ==, sum);

  /* Interleaved alternating big sign */
  double iabs_data[4] = {
    1,
    1e100,
    1,
    -1e100,
  };

  v = gsl_vector_view_array((double *)iabs_data, 4);
  sum = qdm_vector_sum(&v.vector);
  munit_assert_double(2, ==, sum);

  /* Long sequence */
  gsl_vector *ls0_data = qdm_vector_seq(100, 1000, 1);
  sum = qdm_vector_sum(ls0_data);
  munit_assert_double(495550, ==, sum);
  gsl_vector_free(ls0_data);

  gsl_vector *ls1_data = qdm_vector_seq(10, 100, 0.1);
  sum = qdm_vector_sum(ls1_data);
  munit_assert_double_equal(49555, sum, 10);
  gsl_vector_free(ls1_data);

  gsl_vector *ls2_data = qdm_vector_seq(1, 10, 0.01);
  sum = qdm_vector_sum(ls2_data);
  munit_assert_double_equal(4955.5, sum, 10);
  gsl_vector_free(ls2_data);

  gsl_vector *ls3_data = qdm_vector_seq(0.000000001, 0.00000001, 0.00000000001);
  sum = qdm_vector_sum(ls3_data);
  munit_assert_double_equal(4.9555e-06, sum, 10);
  gsl_vector_free(ls3_data);

  return MUNIT_OK;
}

static
MunitResult
test_qdm_matrix_select_upper_triangle(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  double grid_data[] = {
      1,2,3,
      4,5,6,
      7,8,9,
  };

  gsl_matrix_view grid = gsl_matrix_view_array(grid_data, 3, 3);

  double expected_data[] = {
      1,2,3,
      0,5,6,
      0,0,9,
  };

  gsl_matrix_view expected = gsl_matrix_view_array(expected_data, 3, 3);

  qdm_matrix_select_upper_triangle(&grid.matrix);

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      double e = gsl_matrix_get(&expected.matrix, i, j);
      double r = gsl_matrix_get(&grid.matrix, i, j);

      munit_logf(MUNIT_LOG_DEBUG, "(%zu, %zu)", i, j);
      munit_assert_double(e, ==, r);
    }
  }

  return MUNIT_OK;
}

static
MunitResult
test_qdm_matrix_to_csc_matrix(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  int status = 0;

  int n = 2;

  int P_nnz = 3;
  double P_x[3] = {4.0, 1.0, 2.0};
  int P_i[3] = {0, 0, 1};
  int P_p[3] = {0, 1, 3};

  double input_data[] = {
    4, 1,
    1, 2,
  };

  gsl_matrix_view input = gsl_matrix_view_array(input_data, 2, 2);
  qdm_matrix_select_upper_triangle(&input.matrix);

  munit_log(MUNIT_LOG_DEBUG, "input");
  qdm_matrix_csv_fwrite(stderr, &input.matrix);

  csc *result = NULL;
  status = qdm_matrix_to_csc_matrix(&result, &input.matrix);
  munit_assert_int(status, ==, 0);

  munit_assert_int(P_nnz, ==, result->nzmax);
  munit_assert_int(n, ==, result->n);

  for (size_t i = 0; i < (size_t)(result->n + 1); i++) {
    munit_logf(MUNIT_LOG_DEBUG, "result.p[%zu] = %i", i, result->p[i]);
    munit_assert_int(result->p[i], ==, P_p[i]);
  }

  for (size_t i = 0; i < (size_t)result->nzmax; i++) {
    munit_logf(MUNIT_LOG_DEBUG, "result.i[%zu] = %i", i, result->i[i]);
    munit_assert_int(result->i[i], ==, P_i[i]);
  }

  for (size_t i = 0; i < (size_t)result->nzmax; i++) {
    munit_logf(MUNIT_LOG_DEBUG, "result.x[%zu] = %.17g", i, result->x[i]);
    munit_assert_double(result->x[i], ==, P_x[i]);
  }

  return MUNIT_OK;
}

static
MunitResult
test_qdm_vector_hd5_write(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  int status = 0;

  double v_data[] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9,
  };

  gsl_vector_view v = gsl_vector_view_array(v_data, 9);

  int tf_fd = 0;
  char tf_name[] = "/tmp/qdm-test-hd5-vector-XXXXXX";
  tf_fd = mkstemp(tf_name);
  close(tf_fd);

  hid_t tf_file = qdm_data_create_file(tf_name);
  munit_assert_int(tf_file, >=, 0);
  unlink(tf_name);

  hid_t tf_group = qdm_data_create_group(tf_file, "data");
  munit_assert_int(tf_group, >=, 0);

  status = qdm_vector_hd5_write(tf_group, "v", &v.vector);
  munit_assert_int(status, ==, 0);

  status = H5Gclose(tf_group);
  munit_assert_int(status, ==, 0);

  status = H5Fclose(tf_file);
  munit_assert_int(status, ==, 0);

  return MUNIT_OK;
}

static
MunitResult
test_qdm_matrix_hd5_write(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  int status = 0;

  double m_data[] = {
    1, 2, 3,
    4, 5, 6,
    7, 8, 9,
  };

  gsl_matrix_view m = gsl_matrix_view_array(m_data, 3, 3);

  int tf_fd = 0;
  char tf_name[] = "/tmp/qdm-test-hd5-matrix-XXXXXX";
  tf_fd = mkstemp(tf_name);
  close(tf_fd);

  hid_t tf_file = qdm_data_create_file(tf_name);
  munit_assert_int(tf_file, >=, 0);
  unlink(tf_name);

  hid_t tf_group = qdm_data_create_group(tf_file, "data");
  munit_assert_int(tf_group, >=, 0);

  status = qdm_matrix_hd5_write(tf_group, "m", &m.matrix);
  munit_assert_int(status, ==, 0);

  status = H5Gclose(tf_group);
  munit_assert_int(status, ==, 0);

  status = H5Fclose(tf_file);
  munit_assert_int(status, ==, 0);

  return MUNIT_OK;
}

static char *test_param_qdm_vector_seq_min[] = {
  "1",
  "10",
  "0.1",
  NULL,
};

static char *test_param_qdm_vector_seq_max[] = {
  "10",
  "20",
  "1.1",
  NULL,
};

static char *test_param_qdm_vector_seq_delta[] = {
  "1",
  "0.1",
  "2",
  NULL,
};

static MunitParameterEnum test_params_qdm_vector_seq[] = {
  { "qdm_vector_seq_min"   , test_param_qdm_vector_seq_min   }  ,
  { "qdm_vector_seq_max"   , test_param_qdm_vector_seq_max   }  ,
  { "qdm_vector_seq_delta" , test_param_qdm_vector_seq_delta }  ,
  { NULL                   , NULL                            }  ,
};

MunitTest tests_gsl[] = {
  {
    "/qdm_vector_seq"          , // name
    test_qdm_vector_seq        , // test
    NULL                       , // setup
    NULL                       , // tear_down
    MUNIT_TEST_OPTION_NONE     , // options
    test_params_qdm_vector_seq , // parameters
  },
  {
    "/qdm_matrix_det_tmm"   , // name
    test_qdm_matrix_det_tmm , // test
    NULL                    , // setup
    NULL                    , // tear_down
    MUNIT_TEST_OPTION_NONE  , // options
    NULL                    , // parameters
  },
  {
    "/qdm_vector_quantile"   , // name
    test_qdm_vector_quantile , // test
    NULL                     , // setup
    NULL                     , // tear_down
    MUNIT_TEST_OPTION_NONE   , // options
    NULL                     , // parameters
  },
  {
    "/qdm_vector_sum"      , // name
    test_qdm_vector_sum    , // test
    NULL                   , // setup
    NULL                   , // tear_down
    MUNIT_TEST_OPTION_NONE , // options
    NULL                   , // parameters
  },
  {
    "/qdm_matrix_select_upper_triangle"   , // name
    test_qdm_matrix_select_upper_triangle , // test
    NULL                                  , // setup
    NULL                                  , // tear_down
    MUNIT_TEST_OPTION_NONE                , // options
    NULL                                  , // parameters
  },
  {
    "/qdm_matrix_to_csc_matrix"   , // name
    test_qdm_matrix_to_csc_matrix , // test
    NULL                          , // setup
    NULL                          , // tear_down
    MUNIT_TEST_OPTION_NONE        , // options
    NULL                          , // parameters
  },
  {
    "/qdm_vector_hd5_write"   , // name
    test_qdm_vector_hd5_write , // test
    NULL                      , // setup
    NULL                      , // tear_down
    MUNIT_TEST_OPTION_NONE    , // options
    NULL                      , // parameters
  },
  {
    "/qdm_matrix_hd5_write"   , // name
    test_qdm_matrix_hd5_write , // test
    NULL                      , // setup
    NULL                      , // tear_down
    MUNIT_TEST_OPTION_NONE    , // options
    NULL                      , // parameters
  },
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL },
};

static const MunitSuite test_suite_gsl = {
  "/gsl"                  , // name
  tests_gsl               , // tests
  NULL                    , // suites
  1                       , // iterations
  MUNIT_SUITE_OPTION_NONE , // options
};
