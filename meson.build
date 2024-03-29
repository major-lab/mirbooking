project('mirbooking', 'c',
        version: '2.5.0-dev',
        meson_version: '>=0.46.0')

api_version = '2.5'

glib_dep = dependency('glib-2.0', version: '>=2.44')
gobject_dep = dependency('gobject-2.0', version: '>=2.44')
gio_dep = dependency('gio-2.0', version: '>=2.44')
gio_unix_dep = dependency('gio-unix-2.0', version: '>=2.44')
m_dep = meson.get_compiler('c').find_library('m')

if get_option('with_openmp')
    openmp_dep = dependency('openmp')
    add_project_arguments('-DHAVE_OPENMP', language: 'c')
else
    openmp_dep = []
    add_project_arguments('-Wno-unknown-pragmas', language: 'c')
    warning('OpenMP (-Dwith_openmp=true) is recommended for parallelizing the workload.')
endif

if get_option('with_mpi')
    mpi_dep = dependency('mpi', language: 'c')
    add_project_arguments('-DHAVE_MPI', language: 'c')
else
    mpi_dep = []
endif

# Intel MKL
if get_option('with_mkl')
    mkl_root = get_option('mkl_root')
    mkl_arch = get_option('mkl_arch')
    mkl_libdir = join_paths(mkl_root, 'lib', mkl_arch)
    if get_option('with_mkl_ilp64')
        mkl_interface = declare_dependency(compile_args: '-DMKL_ILP64',
                                       dependencies: meson.get_compiler('c').find_library('mkl_intel_ilp64', dirs: mkl_libdir))
    else
        mkl_interface = meson.get_compiler('c').find_library('mkl_intel_lp64', dirs: mkl_libdir)
        warning('MKL 64-bit interface (-Dwith_mkl_ilp64=true) is required for solving large systems.')
    endif
    mkl_core = meson.get_compiler('c').find_library('mkl_core', dirs: mkl_libdir)

    if get_option('with_mkl_tbb')
        tbb = meson.get_compiler('c').find_library('tbb')
        mkl_tbb_thread = meson.get_compiler('c').find_library('mkl_tbb_thread', dirs: mkl_libdir)
        mkl_thread = declare_dependency(dependencies: [mkl_tbb_thread, tbb])
    elif get_option('with_openmp')
        mkl_gnu_thread = meson.get_compiler('c').find_library('mkl_gnu_thread', dirs: mkl_libdir)
        mkl_thread = declare_dependency(dependencies: [mkl_gnu_thread, openmp_dep])
    else
        mkl_thread = meson.get_compiler('c').find_library('mkl_sequential', dirs: mkl_libdir)
        warning('OpenMP (-Dwith_openmp=true) or MKL TBB (-Dwith_mkl_tbb=true) is recommended for parallelizing MKL solvers.')
    endif

    if get_option('with_mpi')
        if get_option('with_mkl_ilp64')
            mkl_blacs = meson.get_compiler('c').find_library('mkl_blacs_openmpi_ilp64', dirs: mkl_libdir)
        else
            mkl_blacs = meson.get_compiler('c').find_library('mkl_blacs_openmpi_lp64', dirs: mkl_libdir)
        endif
    else
        mkl_blacs = []
    endif

    mkl_pthread = meson.get_compiler('c').find_library('pthread')
    mkl_m = meson.get_compiler('c').find_library('m')
    mkl_dl = meson.get_compiler('c').find_library('dl')
    mkl_dep = declare_dependency(compile_args: ['-m64', '-I@0@'.format(join_paths(mkl_root, 'include'))],
                                 link_args: ['-Wl,--no-as-needed'],
                                 dependencies: [mkl_interface, mkl_thread, mkl_core, mkl_blacs, mpi_dep, mkl_pthread, mkl_m, mkl_dl])

    if get_option('with_pardiso')
        error('MKL conflicts with PARDISO solver, use -Dwith_pardiso=false instead.')
    endif

    if get_option('with_mkl_ilp64')
        if not (get_option('with_mkl_cblas') or get_option('with_openblas64'))
            error('MKL 64-bit interface conflicts with 32-bit BLAS, use -Dwith_mkl_cblas=true or -Dwith_openblas64=true instead.')
        endif
        if get_option('with_superlu') or get_option('with_superlu_mt') or get_option('with_umfpack') or get_option('with_pardiso')
            error('MKL 64-bit interface conflicts with the following solvers: LAPACK, SuperLU, SuperLU-MT, UMFPACK and PARDISO.')
        endif
    endif

    if get_option('with_mkl_cluster') and not get_option('with_mpi')
        error('MKL Cluster requires MPI (-Dwith_mpi=true).')
    endif
endif

if not ((get_option('with_mkl') and get_option('with_mkl_dss')) or get_option('with_pardiso') or get_option('with_umfpack'))
    warning('Any of MKL DSS (-Dwith_mkl=true -Dwith_mkl_dss=true), PARDISO (-Dwith_pardiso=true) or UMFPACK (-Dwith_umfpack=true) is recommended for solving sparse systems efficiently.')
endif

if get_option('with_mkl') and get_option('with_mkl_cblas')
    add_project_arguments('-DHAVE_MKL_CBLAS', language: 'c')
    cblas_dep = mkl_dep
elif get_option('with_openblas') or get_option('with_openblas64')
    add_project_arguments('-DHAVE_OPENBLAS', language: 'c')
    if get_option('with_openmp')
        openblas_flavour = 'o'
    else
        openblas_flavour = ''
    endif
    if get_option('with_openblas64')
        openblas_interface = '64'
    else
        openblas_interface = ''
    endif
    cblas_dep = meson.get_compiler('c').find_library('openblas' + openblas_flavour + openblas_interface)
elif get_option('with_atlas')
    cblas_dep = dependency('atlas', required: false)
    if not cblas_dep.found()
        cblas_dep = dependency('blas-atlas')
    endif
else
    cblas_dep = dependency('cblas', required: false)
    if not cblas_dep.found()
        cblas_dep = dependency('blas-netlib')
    endif
endif

if get_option('with_mkl') and get_option('with_mkl_fftw3')
    add_project_arguments('-DHAVE_FFTW3', language: 'c')
    fftw3_dep = declare_dependency(include_directories: include_directories(join_paths(mkl_root, 'include/fftw')),
                                   dependencies: mkl_dep)
elif get_option('with_fftw3')
    add_project_arguments('-DHAVE_FFTW3', language: 'c')
    fftw3_dep = dependency('fftw3')
else
    fftw3_dep = []
    warning('Any of FFTW3 (-Dwith_fftw3=true) or MKL FFTW3 (-Dwith_mkl=true -Dwith_mkl_fftw3=true) is recommended to compute accurate target silencing.')
endif

if get_option('with_mkl')
    lapack_dep = declare_dependency(dependencies: mkl_dep,
                                    compile_args: '-DHAVE_MKL_LAPACK')
else
    lapack_dep = dependency('lapacke', required: false)
    if not lapack_dep.found()
        lapack_dep = meson.get_compiler('c').find_library('lapacke')
    endif
endif


pymod = import('python')
python3 = pymod.find_installation('python3')

subdir('contrib/poisson-binomial')
subdir('contrib/sparse')
subdir('contrib/odeint')
subdir('src')
subdir('bin')
subdir('tools')
subdir('data')
subdir('tests')
