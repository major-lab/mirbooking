project('mirbooking', 'c',
        version: '1.3.0-dev')

glib_dep = dependency('glib-2.0', version: '>=2.44')
gobject_dep = dependency('gobject-2.0', version: '>=2.44')
gio_dep = dependency('gio-2.0', version: '>=2.44')
m_dep = meson.get_compiler('c').find_library('m')

mcff_bin = find_program('mcff')
if mcff_bin.found()
    add_project_arguments(['-DHAVE_MCFF', '-DMCFF="@0@"'.format(mcff_bin.path())], language: 'c')
endif

# Intel MKL
if get_option('with_mkl')
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
else
    warning('MKL is recommended for solving sparse systems efficiently.')
endif

if get_option('with_mkl') and get_option('with_fftw3')
    add_project_arguments('-DHAVE_FFTW3', language: 'c')
    fftw3_dep = declare_dependency(include_directories: include_directories(join_paths(mkl_root, 'include/fftw')),
                                   dependencies: mkl_dep)
elif get_option('with_fftw3')
    add_project_arguments('-DHAVE_FFTW3', language: 'c')
    fftw3_dep = dependency('fftw3', required: get_option('with_fftw3'))
else
    warning('FFTW3 is recommended to compute accurate target silencing.')
endif

if get_option('with_openmp')
    add_project_arguments('-DHAVE_OPENMP', language: 'c')
    openmp_dep = declare_dependency(compile_args: '-fopenmp',
                                    dependencies: meson.get_compiler('c').find_library('gomp'))
else
    add_project_arguments('-Wno-unknown-pragmas', language: 'c')
    openmp_dep = []
    warning('OpenMP is recommended to perform computation in reasonable time.')
endif

if get_option('with_mpi')
    add_project_arguments('-DHAVE_MPI', language: 'c')
    mpicc_bin = find_program('mpicc')
    mpi_dep = declare_dependency(compile_args: run_command(mpicc_bin, ['-showme:compile:']).stdout().split(),
                                 link_args: run_command(mpicc_bin, ['-showme:link:']).stdout().split())
else
    mpi_dep = []
endif

if get_option('with_introspection')
    python3_bin = import('python3').find_python()
    python3_module_deps = ['gi', 'Bio', 'pandas']
    foreach module_dep : python3_module_deps
        result = run_command([python3_bin, '-c', 'import @0@'.format(module_dep)])
        if result.returncode() != 0
            error('Python 3 module \'@0@\' not found: @1@'.format(module_dep, result.stderr()))
        endif
    endforeach
endif

subdir('contrib/poisson-binomial')
subdir('contrib/sparse')
subdir('contrib/odeint')
subdir('src')
subdir('bin')
subdir('data')
subdir('tests')
