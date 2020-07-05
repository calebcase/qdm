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

run_target('cppcheck', command : [
  'cppcheck',
  '--project=' + join_paths(meson.build_root(), 'compile_commands.json'),
  '--enable=warning,performance,portability,information',
  '--error-exitcode=1',
  '--suppress=*:/src/github.com/calebcase/qdm/test/test.c',
  '--suppress=*:/src/github.com/calebcase/qdm/subprojects/munit/munit.c',
  '--suppress=*:/src/github.com/calebcase/qdm/subprojects/munit/munit.h',
])
