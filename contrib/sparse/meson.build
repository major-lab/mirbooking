sparse_solvers_compile_args = []
sparse_solvers_sources = []
sparse_solvers_deps = []

if get_option('with_superlu')
    superlu_dep = meson.get_compiler('c').find_library('superlu')
    sparse_solvers_compile_args += ['-DHAVE_SUPERLU']
    sparse_solvers_sources += ['sparse-superlu.c']
    sparse_solvers_deps += superlu_dep
endif

if get_option('with_superlu_mt')
    superlu_mt_dep = meson.get_compiler('c').find_library('superlumt64_d')
    sparse_solvers_compile_args += ['-DHAVE_SUPERLU_MT', '-D_LONG_INT']
    sparse_solvers_sources += ['sparse-superlu-mt.c']
    sparse_solvers_deps += superlu_mt_dep
endif

if get_option('with_umfpack')
    umfpack_dep = meson.get_compiler('c').find_library('umfpack')
    sparse_solvers_compile_args += ['-DHAVE_UMFPACK']
    sparse_solvers_sources += 'sparse-umfpack.c'
    sparse_solvers_deps += umfpack_dep
endif

if get_option('with_mkl_dss') or get_option('with_mkl_cluster')
    mkl_root = get_option('mkl_root')
    mkl_arch = get_option('mkl_arch')
    mkl_libdir = join_paths(mkl_root, 'lib', mkl_arch)
    mkl_ilp64 = meson.get_compiler('c').find_library('mkl_intel_ilp64', dirs: mkl_libdir)
    mkl_core = meson.get_compiler('c').find_library('mkl_core', dirs: mkl_libdir)
    mkl_gnu_thread = meson.get_compiler('c').find_library('mkl_gnu_thread', dirs: mkl_libdir)
    mkl_gomp = meson.get_compiler('c').find_library('gomp')
    mkl_pthread = meson.get_compiler('c').find_library('pthread')
    mkl_m = meson.get_compiler('c').find_library('m')
    mkl_dl = meson.get_compiler('c').find_library('dl')
    mkl_dep = declare_dependency(compile_args: ['-DMKL_ILP64', '-m64', '-I@0@'.format(join_paths(mkl_root, 'include'))],
                                 link_args: ['-Wl,--no-as-needed'],
                                 dependencies: [mkl_ilp64, mkl_gnu_thread, mkl_core, mkl_gomp, mkl_pthread, mkl_m, mkl_dl])

    if get_option('with_mkl_dss')
        sparse_solvers_compile_args += ['-DHAVE_MKL_DSS']
        sparse_solvers_sources += 'sparse-mkl-dss.c'
        sparse_solvers_deps += mkl_dep
    endif

    if get_option('with_mkl_cluster')
        sparse_solvers_compile_args += ['-DHAVE_MKL_CLUSTER']
        sparse_solvers_sources += 'sparse-mkl-cluster.c'
        mkl_blacs_openmpi_ilp64 = meson.get_compiler('c').find_library('mkl_blacs_openmpi_ilp64', dirs: mkl_libdir)
        sparse_solvers_deps += mkl_blacs_openmpi_ilp64
        sparse_solvers_deps += mpi_dep
    endif
endif

sparse_lib = library('sparse', ['sparse.c'] + sparse_solvers_sources,
                     c_args: sparse_solvers_compile_args,
                     dependencies: sparse_solvers_deps,
                     install: true)
sparse_dep = declare_dependency(include_directories: include_directories('.'),
                                link_with: sparse_lib)

test('sparse', executable('sparse-test', 'sparse-test.c',
                          dependencies: [sparse_dep]))

if get_option('with_mkl_cluster') and get_option('with_mpi')
test('sparse-mpi', executable('sparse-mpi-test', 'sparse-mpi-test.c',
     dependencies: [sparse_dep, mpi_dep]))
endif
