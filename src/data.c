/* Functions for working with the HDF5 files. */

#define _GNU_SOURCE // Required for asprintf.
#include <stdio.h>
#undef _GNU_SOURCE

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
    hid_t id,
    size_t iteration,
    const qdm_intermediate_result *result
)
{
  int status = 0;

  char *group_path = NULL;
  hid_t group = -1;

  status = asprintf(&group_path, "output/intermediate/iteration_%zu", iteration);
  if (status < 0) {
    goto cleanup;
  }

  group = qdm_data_create_group(id, group_path);
  if (group < 0) {
    status = group;
    goto cleanup;
  }

  status = qdm_matrix_hd5_write(group, "theta", result->theta);
  if (status < 0) {
    goto cleanup;
  }

  status = qdm_matrix_hd5_write(group, "theta_star", result->theta_star);
  if (status < 0) {
    goto cleanup;
  }

  status = qdm_vector_hd5_write(group, "ll", result->ll);
  if (status < 0) {
    goto cleanup;
  }

  status = qdm_vector_hd5_write(group, "tau", result->tau);
  if (status < 0) {
    goto cleanup;
  }

  status = qdm_vector_hd5_write(group, "xi", result->xi);
  if (status < 0) {
    goto cleanup;
  }

cleanup:
  if (group >= 0) {
    H5Gclose(group);
  }

  if (group_path != NULL) {
    free(group_path);
  }

  return status;
}

hid_t
qdm_data_create_file(
    const char *path
) {
  return H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

hid_t
qdm_data_create_group(
    hid_t id,
    const char *name
) {
  int status = 0;

  hid_t gcpl  = -1;
  hid_t group = -1;

  gcpl = H5Pcreate(H5P_LINK_CREATE);
  if (gcpl < 0) {
    status = gcpl;
    goto cleanup;
  }

  status = H5Pset_create_intermediate_group(gcpl, 1);
  if (status < 0) {
    goto cleanup;
  }

  group = H5Gcreate(id, name, gcpl, H5P_DEFAULT, H5P_DEFAULT);
  if (group < 0) {
    status = group;
    goto cleanup;
  }

cleanup:
  if (gcpl >= 0) {
    H5Pclose(gcpl);
  }

  if (status == 0) {
    return group;
  }

  return status;
}
