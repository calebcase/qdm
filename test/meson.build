munit_exe = executable('munit', 'test.c',
  include_directories : qdm_inc,
  link_with           : qdm_lib,
  dependencies: [
    dependency('munit', fallback: ['munit', 'munit_dep']),
    dependency('gsl'),
  ],
)

test('munit', munit_exe, args: ['--show-stderr'])
