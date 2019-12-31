#ifndef QDM_TEST_H
#define QDM_TEST_H 1

#define munit_log_size(l, x)           \
  do {                                 \
    munit_logf((l), #x " = %zu", (x)); \
  } while (0)

#define munit_log_vector(l, v)                          \
  do {                                                  \
    munit_log((l), #v);                                 \
    for (size_t i = 0; i < (v)->size; i++) {            \
      munit_logf((l), "%.17g", gsl_vector_get((v), i)); \
    }                                                   \
  } while (0)

#define munit_log_vector_dim(l, v)           \
  do {                                       \
    munit_logf((l), #v " (%zu)", (v)->size); \
  } while (0)

#define munit_log_matrix_dim(l, m)                             \
  do {                                                         \
    munit_logf((l), #m " (%zu, %zu)", (m)->size1, (m)->size2); \
  } while (0)

#define munit_vector_equal(a, b, prec)                   \
  do {                                                   \
    munit_assert_size((a)->size, ==, (b)->size);         \
                                                         \
    for (size_t i = 0; i < (a)->size1; i++) {            \
      double e = gsl_vector_get((a), i, j);              \
      double r = gsl_vector_get((b), i, j);              \
                                                         \
      munit_logf(MUNIT_LOG_DEBUG, "(%zu, %zu)", i, j);   \
      munit_assert_double_equal(e, r, prec);             \
    }                                                    \
  } while (0)

#define munit_matrix_equal(a, b, prec)                   \
  do {                                                   \
    munit_assert_size((a)->size1, ==, (b)->size1);       \
    munit_assert_size((a)->size2, ==, (b)->size2);       \
                                                         \
    for (size_t i = 0; i < (a)->size1; i++) {            \
      for (size_t j = 0; j < (a)->size2; j++) {          \
        double e = gsl_matrix_get((a), i, j);            \
        double r = gsl_matrix_get((b), i, j);            \
                                                         \
        munit_logf(MUNIT_LOG_DEBUG, "(%zu, %zu)", i, j); \
        munit_assert_double_equal(e, r, prec);           \
      }                                                  \
    }                                                    \
  } while (0)

#endif /* QDM_TEST_H */
