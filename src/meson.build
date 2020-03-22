c_args = [
    '-DHAVE_INLINE',
    '-g3',
]

# libqdm

libqdm_src = [
  'data.c',
  'gsl.c',
  'ispline.c',
  'knots.c',
  'logl.c',
  'mspline.c',
  'state.c',
  'theta.c',
]

libqdm_dep = [
  dependency('gsl' , version: '>=2.6' , static: true),
  dependency('hdf5', version: '>=1.10', static: true , language: 'c'),
  dependency('osqp', version: '>=0.6' , static: true),
]

libqdm = static_library('qdm', libqdm_src,
  c_args              : c_args,
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
  dependency('argtable2', version: '>=2.13', static: true),
  dependency('gsl'      , version: '>=2.6' , static: true),
  dependency('hdf5'     , version: '>=1.10', static: true , language: 'c'),
  dependency('inih'     , version: '>=0.47', static: true),
  dependency('osqp'     , version: '>=0.6' , static: true),
]

executable('qdm', qdm_src,
  c_args              : c_args,
  include_directories : libqdm_inc,
  link_with           : libqdm,
  dependencies        : qdm_dep,
  install             : true,

  # Several of these are manually added so that static compilation will work.
  # Most of them are dependencies of hdf5 that are not handled properly with
  # pkg-config. There is a warning during compilation about dlopen, but as far
  # as I can tell the code path for dlopen is never called.
  link_args           : ['-static', '-pthread', '-laec', '-lsz', '-lz', '-ldl'],
)
