/* Functions for working with the HDF5 data file. */
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
