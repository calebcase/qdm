project('qdm', 'c',
  version : '0.1.0',
  license : 'Apache2',
  default_options: [
    'b_coverate=true',
    'warning_level=3',
    'werror=true',
  ],
)

libqdm_inc = include_directories('include')

subdir('include')
subdir('src')
subdir('test')

libqdm_pkg = import('pkgconfig')
libqdm_pkg.generate(libqdm,
  description : 'qdm',
)

qdm_dep = declare_dependency(
  include_directories : libqdm_inc,
  link_with           : libqdm,
)
