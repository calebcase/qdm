#include <munit.h>

#include <qdm.h>

#include "test.h"

static
MunitResult
test_qdm_ijk_get_k(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  double data[] = {
    0, 1, 2, 3,
    1, 1, 2, 3,
    2, 2, 2, 3,

    1, 2, 0, 3,
    2, 2, 0, 3,
    0, 0, 0, 3,

    2, 0, 1, 3,
    0, 0, 1, 3,
    1, 1, 1, 3,

    3, 3, 3, 3,
    3, 3, 3, 3,
    3, 3, 3, 3,
  };

  qdm_ijk_view t = qdm_ijk_view_array(data, 3, 4, 4);

  /* Row 0 */

  {
    double expect_data[] = {
      0, 1, 2, 3,
    };
    gsl_vector_view expect = gsl_vector_view_array(expect_data, 4);

    gsl_vector_view k = qdm_ijk_get_k(&t.ijk, 0, 0);

    munit_vector_equal(&expect.vector, &k.vector, 3);
  }

  {
    double expect_data[] = {
      1, 2, 0, 3,
    };
    gsl_vector_view expect = gsl_vector_view_array(expect_data, 4);

    gsl_vector_view k = qdm_ijk_get_k(&t.ijk, 0, 1);

    munit_vector_equal(&expect.vector, &k.vector, 3);
  }

  {
    double expect_data[] = {
      2, 0, 1, 3,
    };
    gsl_vector_view expect = gsl_vector_view_array(expect_data, 4);

    gsl_vector_view k = qdm_ijk_get_k(&t.ijk, 0, 2);

    munit_vector_equal(&expect.vector, &k.vector, 3);
  }

  {
    double expect_data[] = {
      3, 3, 3, 3,
    };
    gsl_vector_view expect = gsl_vector_view_array(expect_data, 4);

    gsl_vector_view k = qdm_ijk_get_k(&t.ijk, 0, 3);

    munit_vector_equal(&expect.vector, &k.vector, 3);
  }

  /* Row 1 */

  {
    double expect_data[] = {
      1, 2, 0, 3,
    };
    gsl_vector_view expect = gsl_vector_view_array(expect_data, 4);

    gsl_vector_view k = qdm_ijk_get_k(&t.ijk, 1, 0);

    munit_vector_equal(&expect.vector, &k.vector, 3);
  }

  {
    double expect_data[] = {
      1, 2, 0, 3,
    };
    gsl_vector_view expect = gsl_vector_view_array(expect_data, 4);

    gsl_vector_view k = qdm_ijk_get_k(&t.ijk, 1, 1);

    munit_vector_equal(&expect.vector, &k.vector, 3);
  }

  {
    double expect_data[] = {
      2, 0, 1, 3,
    };
    gsl_vector_view expect = gsl_vector_view_array(expect_data, 4);

    gsl_vector_view k = qdm_ijk_get_k(&t.ijk, 1, 2);

    munit_vector_equal(&expect.vector, &k.vector, 3);
  }

  {
    double expect_data[] = {
      3, 3, 3, 3,
    };
    gsl_vector_view expect = gsl_vector_view_array(expect_data, 4);

    gsl_vector_view k = qdm_ijk_get_k(&t.ijk, 1, 3);

    munit_vector_equal(&expect.vector, &k.vector, 3);
  }


  return MUNIT_OK;
}

static
MunitResult
test_qdm_ijk_get_ij(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  double data[] = {
    0, 1, 2, 3,
    1, 1, 2, 3,
    2, 2, 2, 3,

    1, 2, 0, 3,
    2, 2, 0, 3,
    0, 0, 0, 3,

    2, 0, 1, 3,
    0, 0, 1, 3,
    1, 1, 1, 3,

    3, 3, 3, 3,
    3, 3, 3, 3,
    3, 3, 3, 3,
  };

  qdm_ijk_view t = qdm_ijk_view_array(data, 3, 4, 4);

  {
    double expect_data[] = {
      0, 1, 2, 3,
      1, 1, 2, 3,
      2, 2, 2, 3,
    };
    gsl_matrix_view expect = gsl_matrix_view_array(expect_data, 3, 4);

    gsl_matrix_view ij = qdm_ijk_get_ij(&t.ijk, 0);

    munit_matrix_equal(&expect.matrix, &ij.matrix, 3);
  }

  {
    double expect_data[] = {
      1, 2, 0, 3,
      2, 2, 0, 3,
      0, 0, 0, 3,
    };
    gsl_matrix_view expect = gsl_matrix_view_array(expect_data, 3, 4);

    gsl_matrix_view ij = qdm_ijk_get_ij(&t.ijk, 1);

    munit_matrix_equal(&expect.matrix, &ij.matrix, 3);
  }

  {
    double expect_data[] = {
      2, 0, 1, 3,
      0, 0, 1, 3,
      1, 1, 1, 3,
    };
    gsl_matrix_view expect = gsl_matrix_view_array(expect_data, 3, 4);

    gsl_matrix_view ij = qdm_ijk_get_ij(&t.ijk, 2);

    munit_matrix_equal(&expect.matrix, &ij.matrix, 3);
  }

  {
    double expect_data[] = {
      3, 3, 3, 3,
      3, 3, 3, 3,
      3, 3, 3, 3,
    };
    gsl_matrix_view expect = gsl_matrix_view_array(expect_data, 3, 4);

    gsl_matrix_view ij = qdm_ijk_get_ij(&t.ijk, 3);

    munit_matrix_equal(&expect.matrix, &ij.matrix, 3);
  }

  return MUNIT_OK;
}

MunitTest tests_ijk[] = {
  {
    "/qdm_ijk_get_k",       // name
    test_qdm_ijk_get_k,     // test
    NULL,                   // setup
    NULL,                   // tear_down
    MUNIT_TEST_OPTION_NONE, // options
    NULL,                   // parameters
  },
  {
    "/qdm_ijk_get_ij",      // name
    test_qdm_ijk_get_ij,    // test
    NULL,                   // setup
    NULL,                   // tear_down
    MUNIT_TEST_OPTION_NONE, // options
    NULL,                   // parameters
  },
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL },
};

static const MunitSuite test_suite_ijk = {
  "/ijk",                  // name
  tests_ijk,               // tests
  NULL,                    // suites
  1,                       // iterations
  MUNIT_SUITE_OPTION_NONE, // options
};
