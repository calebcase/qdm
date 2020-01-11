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
  dependency('gsl' , version: '>=2.6'),
  dependency('osqp'), # version: '>=0.6'; but osqp doesn't generate a pkg config.
  dependency('hdf5', version: '>=1.10', language: 'c'),
]

libqdm = both_libraries('qdm', libqdm_src,
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
  dependency('argtable2', version: '>=2.13'),
  dependency('inih'     , version: '>=0.47'),
]

executable('qdm', qdm_src,
  c_args              : c_args,
  include_directories : libqdm_inc,
  link_with           : libqdm.get_shared_lib(),
  dependencies        : libqdm_dep + qdm_dep,
)
