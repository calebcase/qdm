munit_exe = executable('munit', 'test.c',
  include_directories : libqdm_inc,
  link_with           : libqdm,
  dependencies: [
    dependency('munit', fallback: ['munit', 'munit_dep']),
    dependency('gsl'),
    dependency('osqp'),
  ],
)

test('munit', munit_exe,
  env: [
    'GSL_RNG_SEED=100',
  ],
  args: [
    '--log-visible', 'debug',
    '--show-stderr',
    '--param', 'y_path', meson.current_source_dir() + '/data/y.csv',
    #'/logl/qdm_logl/within'
    #'/knots/qdm_knots_optimize'
    #'/gsl/qdm_matrix_to_csc_matrix'
  ],
)
