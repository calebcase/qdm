qdm_sources = [
  'gsl.c',
  'ispline.c',
  'knots.c',
  'mspline.c',
]

qdm_dep = [
  dependency('gsl'),
  dependency('argtable2'),
  dependency('jansson'),
]

qdm_lib = both_libraries('qdm', qdm_sources,
  c_args: [
    '-DHAVE_INLINE',
  ],
  include_directories : qdm_inc,
  dependencies        : qdm_dep,
  install             : true,
)

#qdm_dep = [
#  dependency('gsl', static: true),
#  dependency('argtable2', static: true),
#  dependency('jansson', static: true),
#]

#executable('qdm', 'main.c',
#  c_args: [
#    '-DHAVE_INLINE',
#  ],
#  link_args: [
#    '-static',
#    '-g',
#  ],
#  dependencies: [
#    dependency('gsl', static: true),
#    dependency('argtable2', static: true),
#    dependency('jansson', static: true),
#  ],
#)
