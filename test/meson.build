munit_exe = executable('munit', 'test.c',
  include_directories : libqdm_inc,
  link_with           : libqdm,
  dependencies: [
    dependency('gsl',      version: '>=2.6'),
    dependency('hdf5',     version: '>=1.10',   language: 'c'),
    dependency('openblas', version: '>=0.3.7'),
    dependency('osqp',     version: '>=0.6'),

    dependency('munit', fallback: ['munit', 'munit_dep']),
  ],
)

test('munit', munit_exe,
  env: [
    'GSL_RNG_SEED=100',
  ],
  args: [
    '--color', 'always',
  ],
)
