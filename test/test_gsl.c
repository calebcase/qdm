#include <munit.h>
#include <qdm.h>

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

    munit_logf(MUNIT_LOG_DEBUG, "%f", found);
    munit_assert_double_equal(expecting, found, 3);
  }

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
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL },
};

static const MunitSuite test_suite_gsl = {
  "/gsl"                  , // name
  tests_gsl               , // tests
  NULL                    , // suites
  1                       , // iterations
  MUNIT_SUITE_OPTION_NONE , // options
};
