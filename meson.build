project('qdm', 'c',
  version : '0.1.0',
  license : 'Apache2',
  default_options: [
    'b_coverate=true',
    'c_std=c18',
    'warning_level=3',
    'werror=true',
  ],
)

qdm_inc = include_directories('include')

subdir('include')
subdir('src')
subdir('test')

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  libraries   : '0.1',
  name        : 'libqdm',
  filebase    : 'qdm',
  description : 'qdm',
)
