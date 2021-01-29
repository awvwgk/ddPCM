project(
  'ddcosmo',
  'fortran',
  meson_version: '>=0.53',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# ddCOSMO library target
ddcosmo_lib = library(
  meson.project_name(),
  sources: srcs,
  install: install,
)

# Export dependency for other projects and test suite
ddcosmo_inc = ddcosmo_lib.private_dir_include()
ddcosmo_dep = declare_dependency(
  link_with: ddcosmo_lib,
  include_directories: ddcosmo_inc,
)

# Add example driver
subdir('example')

# Package the license files
ddcosmo_lic = files(
  'LICENSE',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    ddcosmo_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name()
  )

  pkg = import('pkgconfig')
  pkg.generate(
    ddcosmo_lib,
    description: 'A fast domain decomposition based implementation of the COSMO solvation model',
  )
endif
