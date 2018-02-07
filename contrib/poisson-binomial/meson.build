pb_deps = [m_dep]
if get_option('with_fftw3')
    pb_deps += fftw3_dep
endif
pb_lib = static_library('pb', 'pb.c', dependencies: pb_deps)
pb_dep = declare_dependency(include_directories: include_directories('.'),
                            link_with: pb_lib)
test('pb', executable('pb-test', 'pb-test.c', dependencies: [pb_dep, glib_dep]))
