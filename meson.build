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

if get_option('with_gomp')
    add_project_arguments('-DHAVE_OPENMP', language: 'c')
    openmp_dep = declare_dependency(compile_args: '-fopenmp',
                                    dependencies: meson.get_compiler('c').find_library('gomp'))
else
    add_project_arguments('-Wno-unknown-pragmas', language: 'c')
    openmp_dep = []
    warning('OpenMP is recommended to perform computation in reasonable time.')
endif

# Intel MPI
if get_option('with_impi')
    add_project_arguments('-DHAVE_MPI', language: 'c')
    impi_root = get_option('impi_root')
    impi_arch = get_option('impi_arch')
    impi_link_args = [
    '-L/soft/bioinfo/linux_RH7/intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/lib/release_mt',
    '-L/soft/bioinfo/linux_RH7/intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/lib',
    '-Xlinker', '--enable-new-dtags', '-Xlinker', '-rpath', '-Xlinker',
    '/soft/bioinfo/linux_RH7/intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/lib/release_mt',
    '-Xlinker', '-rpath', '-Xlinker',
    '/soft/bioinfo/linux_RH7/intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/lib',
    '-Xlinker', '-rpath', '-Xlinker',
    '/opt/intel/mpi-rt/2107.0.0/intel64/lib/release_mt', '-Xlinker', '-rpath',
    '-Xlinker', '/opt/intel/mpi-rt/2017.0.0/intel64/lib', '-lmpifort', '-lmpi',
    '-lmpigi', '-ldl', '-lrt', '-lpthread']

    mpi_dep = declare_dependency(include_directories: include_directories(join_paths(impi_root, impi_arch, 'include')),
                                 link_args: impi_link_args)
elif get_option('with_openmpi')
    add_project_arguments('-DHAVE_MPI', language: 'c')
    mpicc = find_program('mpicc')
    mpi_dep = declare_dependency(compile_args: run_command(mpicc, ['-showme:compile:']).stdout().split(),
                                 link_args: run_command(mpicc, ['-showme:link:']).stdout().split())
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
        warning('MKL 64-bit interface is required for solving large systems.')
    endif
    mkl_core = meson.get_compiler('c').find_library('mkl_core', dirs: mkl_libdir)

    if get_option('with_gomp')
        mkl_thread = meson.get_compiler('c').find_library('mkl_gnu_thread', dirs: mkl_libdir)
    else
        mkl_thread = meson.get_compiler('c').find_library('mkl_sequential', dirs: mkl_libdir)
    endif

    if get_option('with_impi')
        if get_option('with_mkl_ilp64')
            mkl_blacs = meson.get_compiler('c').find_library('mkl_blacs_intelmpi_ilp64', dirs: mkl_libdir)
        else
            mkl_blacs = meson.get_compiler('c').find_library('mkl_blacs_intelmpi_lp64', dirs: mkl_libdir)
        endif
    elif get_option('with_openmpi')
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
                                 dependencies: [mkl_interface, mkl_thread, mkl_core, mkl_blacs, mpi_dep, openmp_dep, mkl_pthread, mkl_m, mkl_dl])
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

python3_bin = import('python3').find_python()
python3_module_deps = ['numpy', 'pandas', 'Bio', 'tqdm']

if get_option('with_introspection')
    python3_module_deps += ['gi']
endif

foreach module_dep : python3_module_deps
    result = run_command([python3_bin, '-c', 'import @0@'.format(module_dep)])
    if result.returncode() != 0
        error('Python 3 module \'@0@\' not found: @1@'.format(module_dep, result.stderr()))
    endif
endforeach

subdir('contrib/poisson-binomial')
subdir('contrib/sparse')
subdir('contrib/odeint')
subdir('src')
subdir('bin')
subdir('data')
subdir('tests')
