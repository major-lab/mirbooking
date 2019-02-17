sparse_solvers_compile_args = []
sparse_solvers_sources = []
sparse_solvers_deps = []

if get_option('with_superlu')
    superlu_dep = meson.get_compiler('c').find_library('superlu')
    if meson.get_compiler('c').has_header('SuperLU/slu_ddefs.h')
    elif meson.get_compiler('c').has_header('superlu/slu_ddefs.h')
        sparse_solvers_compile_args += ['-DHAVE_SUPERLU_LOWERCASE_INCDIR']
    else
        error('Could not find SuperLU header.')
    endif
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

if get_option('with_mkl')
    if get_option('with_mkl_dss')
        sparse_solvers_compile_args += ['-DHAVE_MKL_DSS']
        sparse_solvers_sources += 'sparse-mkl-dss.c'
        sparse_solvers_deps += mkl_dep
    endif

    if get_option('with_mkl_cluster')
        sparse_solvers_compile_args += ['-DHAVE_MKL_CLUSTER']
        sparse_solvers_sources += 'sparse-mkl-cluster.c'
        sparse_solvers_deps += mkl_dep
    endif
endif

if get_option('with_cusolver')
    cuda_root = get_option('cuda_root')
    cudart_lib = meson.get_compiler('c').find_library('cudart', dirs: join_paths(cuda_root, 'lib64'))
    cusparse_lib = meson.get_compiler('c').find_library('cusparse', dirs: join_paths(cuda_root, 'lib64'))
    cusolver_lib = meson.get_compiler('c').find_library('cusolver', dirs: join_paths(cuda_root, 'lib64'))
    cusolver_dep = declare_dependency(include_directories: include_directories(join_paths(cuda_root, 'include')),
                                      dependencies: [cudart_lib, cusparse_lib, cusolver_lib])
    sparse_solvers_compile_args += ['-DHAVE_CUSOLVER']
    sparse_solvers_sources += 'sparse-cusolver.c'
    sparse_solvers_deps += cusolver_dep
endif

if get_option('with_pardiso')
    pardiso_openmp_dep = dependency('openmp')
    gfortran_lib = meson.get_compiler('c').find_library('gfortran')
    pardiso_lib = meson.get_compiler('c').find_library('pardiso')
    pardiso_dep = declare_dependency(dependencies: [gfortran_lib, pardiso_lib, pardiso_openmp_dep])
    sparse_solvers_sources += 'sparse-pardiso.c'
    sparse_solvers_deps += pardiso_dep
    sparse_solvers_compile_args += ['-DHAVE_PARDISO']
endif

sparse_lib = static_library('sparse', ['sparse.c'] + sparse_solvers_sources,
                     c_args: sparse_solvers_compile_args,
                     dependencies: sparse_solvers_deps)
sparse_dep = declare_dependency(include_directories: include_directories('.'),
                                link_with: sparse_lib)

test('sparse', executable('sparse-test', 'sparse-test.c',
                          dependencies: [sparse_dep]))

if get_option('with_mkl') and get_option('with_mkl_cluster')
    test('sparse-mpi', executable('sparse-mpi-test', 'sparse-mpi-test.c',
         dependencies: [sparse_dep, mpi_dep]))
endif
