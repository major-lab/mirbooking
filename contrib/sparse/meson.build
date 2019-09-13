sparse_solvers_compile_args = []
sparse_solvers_modules = []

if get_option('with_superlu')
    superlu_dep = meson.get_compiler('c').find_library('superlu')
    if meson.get_compiler('c').has_header('SuperLU/slu_ddefs.h')
    elif meson.get_compiler('c').has_header('superlu/slu_ddefs.h')
        sparse_solvers_compile_args += ['-DHAVE_SUPERLU_LOWERCASE_INCDIR']
    else
        error('Could not find SuperLU header.')
    endif
    sparse_solvers_compile_args += ['-DHAVE_SUPERLU']
    sparse_superlu_lib = static_library('sparse-superlu', 'sparse-superlu.c',
                                        dependencies: [superlu_dep])
    sparse_solvers_modules += sparse_superlu_lib
endif

if get_option('with_superlu_mt')
    superlu_mt_dep = meson.get_compiler('c').find_library('superlumt64_d')
    sparse_solvers_compile_args += ['-DHAVE_SUPERLU_MT', '-D_LONG_INT']
    sparse_superlu_mt_lib = static_library('sparse-superlu-mt', 'sparse-superlu-mt.c',
                                           dependencies: [superlu_mt_dep])
    sparse_solvers_modules += sparse_superlu_mt_lib
endif

if get_option('with_umfpack')
    umfpack_dep = meson.get_compiler('c').find_library('umfpack')
    sparse_solvers_compile_args += ['-DHAVE_UMFPACK']
    sparse_umfpack_lib = static_library ('sparse-umfpack', 'sparse-umfpack.c',
                                         dependencies: umfpack_dep)
    sparse_solvers_modules += sparse_umfpack_lib
endif

if get_option('with_mkl_dss')
    sparse_solvers_compile_args += ['-DHAVE_MKL_DSS']
    sparse_mkl_dss_lib = static_library('sparse-mkl-dss', 'sparse-mkl-dss.c',
                                        dependencies: mkl_dep)
    sparse_solvers_modules += sparse_mkl_dss_lib
endif

if get_option('with_mkl_cluster')
    sparse_solvers_compile_args += ['-DHAVE_MKL_CLUSTER']
    sparse_mkl_cluster_lib = static_library('sparse-mkl-cluster', 'sparse-mkl-cluster.c',
                                            dependencies: mkl_dep)
    sparse_solvers_modules += sparse_mkl_cluster_lib
endif

if get_option('with_mkl_lapack')
    sparse_solvers_compile_args += ['-DHAVE_MKL_LAPACK']
    sparse_mkl_lapack_lib = static_library('sparse-mkl-lapack', 'sparse-mkl-lapack.c',
        dependencies: mkl_dep)
    sparse_solvers_modules += sparse_mkl_lapack_lib

endif

if get_option('with_cuda') != '' and get_option('with_cusolver')
    cudart_dep = dependency('cudart-@0@'.format(cuda_toolkit_api_version))
    cusparse_dep = dependency('cusparse-@0@'.format(cuda_toolkit_api_version))
    cusolver_dep = dependency('cusolver-@0@'.format(cuda_toolkit_api_version))
    sparse_solvers_compile_args += ['-DHAVE_CUSOLVER']
    sparse_cusolver_lib = static_library('sparse-cusolver', 'sparse-cusolver.c',
                                         dependencies: [cudart_dep, cusparse_dep, cusolver_dep])
    sparse_solvers_modules += sparse_cusolver_lib
endif

if get_option('with_pardiso')
    pardiso_openmp_dep = dependency('openmp')
    gfortran_lib = meson.get_compiler('c').find_library('gfortran')
    pardiso_lib = meson.get_compiler('c').find_library('pardiso')
    sparse_solvers_compile_args += ['-DHAVE_PARDISO']
    sparse_pardiso_lib = static_library('sparse-pardiso', 'sparse-pardiso.c',
        dependencies: [gfortran_lib, pardiso_lib, pardiso_openmp_dep])
    sparse_solvers_modules += sparse_pardiso_lib
endif

sparse_lib = static_library('sparse', 'sparse.c',
                            c_args: sparse_solvers_compile_args,
                            link_with: sparse_solvers_modules)
sparse_dep = declare_dependency(include_directories: include_directories('.'),
                                link_with: sparse_lib)

if get_option('with_superlu')
    test('sparse', executable('sparse-test', 'sparse-test.c',
                              dependencies: sparse_dep))
endif

if get_option('with_mkl') and get_option('with_mkl_cluster')
    test('sparse-mpi', executable('sparse-mpi-test', 'sparse-mpi-test.c',
         dependencies: [sparse_dep, mpi_dep]))
endif
