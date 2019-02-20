project('mirbooking', 'c',
        version: '1.2.4')

glib_dep = dependency('glib-2.0', version: '>=2.44')
gobject_dep = dependency('gobject-2.0', version: '>=2.44')
gio_dep = dependency('gio-2.0', version: '>=2.44')
m_dep = meson.get_compiler('c').find_library('m')

pymod = import('python')
python2 = pymod.find_installation('python2')
python3 = pymod.find_installation('python3')

subdir('src')
subdir('bin')
subdir('data')
subdir('tests')
