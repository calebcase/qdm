#ifndef QDM_DATA_H
#define QDM_DATA_H 1

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

#endif /* QDM_DATA_H */
