# libqdm

libqdm_src = [
  'gsl.c',
  'ispline.c',
  'knots.c',
  'logl.c',
  'mspline.c',
  'state.c',
  'theta.c',
]

libqdm_dep = [
  dependency('gsl', version: '>=2.6'),
  dependency('osqp'),
]

libqdm = both_libraries('qdm', libqdm_src,
  c_args: [
    '-DHAVE_INLINE',
    '-D_POSIX_C_SOURCE=200809L',
    '-g3'
  ],
  include_directories : libqdm_inc,
  dependencies        : libqdm_dep,
  install             : true,
)

# qdm

qdm_src = [
  'config.c',
  'main.c',
]

qdm_dep = [
  dependency('gsl'      , static: true),
  dependency('argtable2', static: true),
  dependency('jansson'  , static: true),
  dependency('inih'     , static: true),
]

executable('qdm', qdm_src,
  c_args: [
    '-DHAVE_INLINE',
    '-g3'
  ],
  link_args: [
    '-static',
    '-g3',
  ],
  dependencies: qdm_dep,
)
