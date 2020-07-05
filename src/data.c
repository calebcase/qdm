/* Functions for working with the HDF5 files. */

#include <hdf5.h>
#include <hdf5_hl.h>

#include <qdm.h>

int
qdm_double_write(
    hid_t id,
    const char *name,
    double value
)
{
  hsize_t dims[1] = {1};
  double data[1] = {value};

  return H5LTmake_dataset_double(id, name, 1, dims, data);
}

int
qdm_int_write(
    hid_t id,
    const char *name,
    int value
)
{
  hsize_t dims[1] = {1};
  int data[1] = {value};

  return H5LTmake_dataset_int(id, name, 1, dims, data);
}

hid_t
qdm_data_create_file(
    const char *path
) {
  hid_t id = -1;

  id = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);
  if (id < 0) {
    id = H5Fcreate(path, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  }

  return id;
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
