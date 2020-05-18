#ifndef QDM_DATA_H
#define QDM_DATA_H 1

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <hdf5.h>

typedef struct {
  double year;
  double month;
  double day;
  double value;
} qdm_data_model;

int
qdm_data_model_read(
    qdm_data_model **records,
    size_t *nrecords,
    const char *file_path,
    const char *group_path,
    const char *table_name
);

typedef int (*qdm_data_model_selector)(const qdm_data_model *record, const void *arg);

size_t
qdm_data_model_selected(
    const qdm_data_model *records,
    size_t nrecords,
    qdm_data_model_selector select,
    const void *arg
);

gsl_vector *
qdm_data_model_year_vector(
    const qdm_data_model *records,
    size_t nrecords,
    qdm_data_model_selector select,
    const void *arg
);

gsl_vector *
qdm_data_model_value_vector(
    const qdm_data_model *records,
    size_t nrecords,
    qdm_data_model_selector select,
    const void *arg
);

int
qdm_data_model_select_by_month(
    const qdm_data_model *record,
    const double *month
);

typedef struct {
  gsl_matrix *theta;
  gsl_matrix *theta_star;

  gsl_vector *ll;
  gsl_vector *tau;

  gsl_vector *xi;
} qdm_intermediate_result;

void
qdm_data_intermediate_free(
    qdm_intermediate_result *result
);

int
qdm_data_intermediate_result_write(
    hid_t id,
    size_t iteration,
    const qdm_intermediate_result *result
);

hid_t
qdm_data_create_file(
    const char *path
);

hid_t
qdm_data_create_group(
    hid_t id,
    const char *name
);

#endif /* QDM_DATA_H */
