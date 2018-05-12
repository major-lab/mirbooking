superlu_dep = meson.get_compiler('c').find_library('superlu')
superlu_mt_dep = meson.get_compiler('c').find_library('superlumt64_d')
sparse_lib = library('sparse', ['sparse.c', 'sparse-superlu.c', 'sparse-superlu-mt.c'], dependencies: [m_dep, superlu_dep, superlu_mt_dep])
sparse_dep = declare_dependency(include_directories: include_directories('.'),
                                link_with: sparse_lib)

test('sparse', executable('sparse-test', 'sparse-test.c',
                          dependencies: [sparse_dep]))
