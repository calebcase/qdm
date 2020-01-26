/* Functions for working with the HDF5 files. */
#include <math.h>
#include <stdlib.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <qdm.h>

#define QDM_DATA_MODEL_NFIELDS (hsize_t)4

int
qdm_data_model_read(
    qdm_data_model **records,
    size_t *nrecords,
    const char *file_path,
    const char *group_path,
    const char *table_name
)
{
  int status = 0;

  hid_t file_id;
  hid_t group_id;

  // FIXME: Check for negative return indicating failure to open.
  file_id = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);
  group_id = H5Gopen(file_id, group_path, H5P_DEFAULT);

  hsize_t nfields = 0;
  hsize_t hnrecords = 0;

  status = H5TBget_table_info(group_id, table_name, &nfields, &hnrecords);
  if (status != 0) {
    goto cleanup;
  }
  *nrecords = hnrecords;

  if (nfields != QDM_DATA_MODEL_NFIELDS) {
    status = 1;

    goto cleanup;
  }

  *records = calloc(*nrecords, sizeof(qdm_data_model));

  const size_t qdm_data_model_offsets[QDM_DATA_MODEL_NFIELDS] = {
    HOFFSET(qdm_data_model, year),
    HOFFSET(qdm_data_model, month),
    HOFFSET(qdm_data_model, day),
    HOFFSET(qdm_data_model, value),
  };

  const size_t qdm_data_model_sizes[QDM_DATA_MODEL_NFIELDS] = {
    sizeof((*records[0]).year),
    sizeof((*records[0]).month),
    sizeof((*records[0]).day),
    sizeof((*records[0]).value),
  };

  status = H5TBread_table(
      group_id,
      table_name,
      sizeof(qdm_data_model),
      qdm_data_model_offsets,
      qdm_data_model_sizes,
      *records
  );
  if (status != 0) {
    free(*records);

    goto cleanup;
  }

cleanup:
  H5Gclose(group_id);
  H5Fclose(file_id);

  return status;
}

size_t
qdm_data_model_selected(
    const qdm_data_model *records,
    size_t nrecords,
    int (*select)(const qdm_data_model *record, void *arg),
    void *arg
)
{
  size_t size = 0;

  for (size_t i = 0; i < nrecords; i++) {
    if (select(&records[i], arg) == 1) {
      size++;
    }
  }

  return size;
}

gsl_vector *
qdm_data_model_value_vector(
    const qdm_data_model *records,
    size_t nrecords,
    int (*select)(const qdm_data_model *record, void *arg),
    void *arg
)
{
  size_t size = qdm_data_model_selected(records, nrecords, select, arg);

  // Allocate the vector and fill it in.
  gsl_vector *v = gsl_vector_alloc(size);

  size_t j = 0;
  for (size_t i = 0; i < nrecords; i++) {
    if (select(&records[i], arg) == 1) {
      gsl_vector_set(v, j, records[i].value);
      j++;
    }
  }

  return v;
}

gsl_vector *
qdm_data_model_year_vector(
    const qdm_data_model *records,
    size_t nrecords,
    int (*select)(const qdm_data_model *record, void *arg),
    void *arg
)
{
  size_t size = qdm_data_model_selected(records, nrecords, select, arg);

  // Allocate the vector and fill it in.
  gsl_vector *v = gsl_vector_alloc(size);

  size_t j = 0;
  for (size_t i = 0; i < nrecords; i++) {
    if (select(&records[i], arg) == 1) {
      gsl_vector_set(v, j, records[i].year);
      j++;
    }
  }

  return v;
}

int
qdm_data_model_select_by_month(
    const qdm_data_model *record,
    double *month
)
{
  if (isnan(record->value)) {
    return 0;
  }

  if (record->month == *month) {
    return 1;
  }

  return 0;
}

int
qdm_data_intermediate_result_write(
    const char *file_path,
    const char *group_path,
    const qdm_intermediate_result *result
)
{
  int status = 0;

  status = qdm_matrix_hd5_write(file_path, group_path, "theta", result->theta);
  if (status < 0) {
    return status;
  }

  status = qdm_matrix_hd5_write(file_path, group_path, "theta_star", result->theta_star);
  if (status < 0) {
    return status;
  }

  status = qdm_vector_hd5_write(file_path, group_path, "ll", result->ll);
  if (status < 0) {
    return status;
  }

  status = qdm_vector_hd5_write(file_path, group_path, "tau", result->tau);
  if (status < 0) {
    return status;
  }

  status = qdm_vector_hd5_write(file_path, group_path, "xi", result->xi);
  if (status < 0) {
    return status;
  }

  return status;
}
