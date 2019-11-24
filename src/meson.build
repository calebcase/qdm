executable('qdm', 'main.c',
  c_args: [
    '-DHAVE_INLINE',
  ],
  link_args: [
    '-static',
    '-g',
  ],
  dependencies: [
    dependency('gsl', static: true),
    dependency('argtable2', static: true),
    dependency('jansson', static: true),
  ],
)
