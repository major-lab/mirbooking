project('miRBooking', 'c',
        version: '1.0.0-dev')

glib_dep = dependency('glib-2.0', version: '>=2.44')
gobject_dep = dependency('gobject-2.0', version: '>=2.44')
gio_dep = dependency('gio-2.0', version: '>=2.44')
m_dep = meson.get_compiler('c').find_library('m')

subdir('src')
subdir('bin')
subdir('tests')
