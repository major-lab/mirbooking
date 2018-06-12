project('mirbooking', 'c',
        version: '1.3.0-dev')

glib_dep = dependency('glib-2.0', version: '>=2.44')
gobject_dep = dependency('gobject-2.0', version: '>=2.44')
gio_dep = dependency('gio-2.0', version: '>=2.44')
m_dep = meson.get_compiler('c').find_library('m')
fftw3_dep = dependency('fftw3', required: get_option('with_fftw3'))
if get_option('with_openmp')
    add_project_arguments('-DHAVE_OPENMP', language: 'c')
    openmp_dep = declare_dependency(compile_args: '-fopenmp',
                                    dependencies: meson.get_compiler('c').find_library('gomp'))
else
    add_project_arguments('-Wno-unknown-pragmas', language: 'c')
    openmp_dep = []
    warning('OpenMP is suggested to perform computation in reasonable time.')
endif
if get_option('with_fftw3')
    add_project_arguments('-DHAVE_FFTW3', language: 'c')
else
    warning('FFTW3 is needed to compute accurate target silencing.')
endif

if get_option('with_mpi')
    add_project_arguments('-DHAVE_MPI', language: 'c')
    mpicc_bin = find_program('mpicc')
    mpi_dep = declare_dependency(compile_args: run_command(mpicc_bin, ['-showme:compile:']).stdout().split(),
                                 link_args: run_command(mpicc_bin, ['-showme:link:']).stdout().split())
else
    mpi_dep = []
    warning('MPI is needed to use the MKL Parallel Direct Sparse Solver for Clusters backend.')
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
subdir('src')
subdir('bin')
subdir('data')
subdir('tests')
