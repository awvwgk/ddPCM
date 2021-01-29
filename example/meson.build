ddcosmo_example_exe = executable(
  '@0@-example'.format(meson.project_name()),
  sources: files('main.f90'),
  dependencies: ddcosmo_dep,
  install: false,
)

foreach t : ['ala25', 'benzene', 'caffeine', 'crambin', 'formaldehyde', 'taxol']
  test(t, ddcosmo_example_exe, args: [files('@0@.txt'.format(t))])
endforeach
