munit_exe = executable('munit', 'test.c',
  include_directories : libqdm_inc,
  link_with           : libqdm,
  dependencies: [
    dependency('gsl'),
    dependency('hdf5', version: '>=1.10', language: 'c'),
    dependency('munit', fallback: ['munit', 'munit_dep']),
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
  ],
)
