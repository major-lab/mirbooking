project('mirbooking', 'c',
        version: '2.3.0-dev',
        default_options: 'c_std=gnu11',
        meson_version: '>=0.46.0')

api_version = '2.3'

glib_dep = dependency('glib-2.0', version: '>=2.44')
gobject_dep = dependency('gobject-2.0', version: '>=2.44')
gio_dep = dependency('gio-2.0', version: '>=2.44')
m_dep = meson.get_compiler('c').find_library('m')

if get_option('with_openmp')
    openmp_dep = dependency('openmp')
    add_project_arguments('-DHAVE_OPENMP', language: 'c')
else
    openmp_dep = []
    add_project_arguments('-Wno-unknown-pragmas', language: 'c')
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
        warning('MKL 64-bit interface is required for solving large systems.')
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
else
    warning('MKL is recommended for solving sparse systems efficiently.')
endif

if get_option('with_mkl')
    add_project_arguments('-DHAVE_MKL_CBLAS', language: 'c')
    cblas_dep = mkl_dep
else
    cblas_dep = dependency('cblas')
endif

if get_option('with_mkl') and get_option('with_mkl_fftw3')
    add_project_arguments('-DHAVE_FFTW3', language: 'c')
    fftw3_dep = declare_dependency(include_directories: include_directories(join_paths(mkl_root, 'include/fftw')),
                                   dependencies: mkl_dep)
    message('MKL FFTW3 will be used.')
elif get_option('with_fftw3')
    add_project_arguments('-DHAVE_FFTW3', language: 'c')
    fftw3_dep = dependency('fftw3')
else
    fftw3_dep = []
    warning('FFTW3 is recommended to compute accurate target silencing.')
endif

pymod = import('python')
python3 = pymod.find_installation('python3')

cuda_toolkit_api_version = get_option('with_cuda')

subdir('contrib/poisson-binomial')
subdir('contrib/sparse')
subdir('contrib/odeint')
subdir('src')
subdir('bin')
subdir('tools')
subdir('data')
subdir('tests')
