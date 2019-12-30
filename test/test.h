#ifndef QDM_TEST_H
#define QDM_TEST_H 1

#define munit_log_size(l, x) \
  do { \
    munit_logf((l), #x " = %zu", (x)); \
  } while (0)

#define munit_log_vector(l, v) \
  do { \
    munit_log((l), #v); \
    for (size_t i = 0; i < (v)->size; i++) { \
      munit_logf((l), "%.17g", gsl_vector_get((v), i)); \
    } \
  } while (0)

#define munit_log_vector_dim(l, v) \
  do { \
    munit_logf((l), #v " (%zu)", (v)->size); \
  } while (0)

#define munit_log_matrix_dim(l, m) \
  do { \
    munit_logf((l), #m " (%zu, %zu)", (m)->size1, (m)->size2); \
  } while (0)

#endif /* QDM_TEST_H */
