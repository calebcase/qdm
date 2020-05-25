#ifndef QDM_DATA_H
#define QDM_DATA_H 1

#include <hdf5.h>

int
qdm_double_write(
    hid_t id,
    const char *name,
    double value
);

int
qdm_int_write(
    hid_t id,
    const char *name,
    int value
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
