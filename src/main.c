#include <stdio.h>

#include <argtable2.h>
#include <qdm.h>

#include "config.h"

int
qdm_main(qdm_config *cfg)
{
  int status = 0;

  qdm_data_model *records = NULL;
  size_t nrecords = 0;

  qdm_config_fprint(stderr, cfg);

  status = qdm_data_model_read(
      &records,
      &nrecords,
      cfg->data_path,
      "input/model",
      "m1"
  );

  for (size_t i = 0; i < nrecords; i++) {
    fprintf(stderr, "record[%zu]: %.0f/%.0f/%.0f %f\n",
        i,
        records[i].year,
        records[i].month,
        records[i].day,
        records[i].value
    );
  }

  return status;
}

int
main(int argc, char **argv)
{
  int status = 0;

  const char *progname = "qdm";
  int nerrors = 0;

  struct arg_file *config_path;
  struct arg_lit *help;
  struct arg_end *end;

  void *argtable[] = {
    config_path = arg_file0("c", "config", "PATH", "path to config"),
    help        = arg_lit0(NULL, "help", "display this help and exit"),
    end         = arg_end(20),
  };

  if (arg_nullcheck(argtable) != 0) {
    fprintf(stderr, "%s: insufficient memory\n",progname);

    status = 1;
    goto cleanup;
  }

  /* Set default values. */
  config_path->filename[0] = "qdm.ini";

  /* Parse command line. */
  nerrors = arg_parse(argc, argv, argtable);

  if (help->count > 0) {
    printf("Usage: %s", progname);
    arg_print_syntax(stdout, argtable, "\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");

    status = 0;
    goto cleanup;
  }

  if (nerrors > 0) {
    arg_print_errors(stderr, end, progname);
    fprintf(stderr, "Try '%s --help' for more information.\n", progname);

    status = 1;
    goto cleanup;
  }

  FILE *config_file = fopen(config_path->filename[0], "r");
  if (config_file == NULL) {
    perror(config_path->filename[0]);

    status = 1;
    goto cleanup;
  }

  qdm_config cfg;
  status = qdm_config_new(config_file, &cfg);
  if (status != 0) {
    fprintf(stderr, "Failed to parse config file.\n");

    status = 1;
    goto cleanup;
  }

  status = qdm_main(&cfg);
  if (status != 0) {
    fprintf(stderr, "Failed to execute qdm routine.\n");

    goto cleanup;
  }

cleanup:
  arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

  return -status;
}
