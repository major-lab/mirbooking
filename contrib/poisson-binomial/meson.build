pb_lib = static_library('pb', 'pb.c', dependencies: [fftw3_dep, m_dep])
pb_dep = declare_dependency(include_directories: include_directories('.'),
                            link_with: pb_lib)
test('pb', executable('pb-test', 'pb-test.c', dependencies: [pb_dep, glib_dep]))
