project('mirbooking', 'c',
        version: '1.1.0')

glib_dep = dependency('glib-2.0', version: '>=2.44')
gobject_dep = dependency('gobject-2.0', version: '>=2.44')
gio_dep = dependency('gio-2.0', version: '>=2.44')
m_dep = meson.get_compiler('c').find_library('m')

python3_bin = import('python3').find_python()
python3_module_deps = [
    'Bio',
    'pandas']
if get_option('with_introspection')
    python3_module_deps += 'gi'
endif
foreach module_dep : python3_module_deps
    if run_command([python3_bin, '-c', 'import @0@'.format(module_dep)]).returncode() != 0
        error('Python 3 module \'@0@\' not found: @1@'.format(module_dep, result.stderr()))
    endif
endforeach

subdir('src')
subdir('bin')
subdir('tests')
