mirbooking_headers = [
    'mirbooking.h',
    'mirbooking-error.h',
    'mirbooking-broker.h',
    'mirbooking-mirna.h',
    'mirbooking-occupant.h',
    'mirbooking-precomputed-score-table.h',
    'mirbooking-score-table.h',
    'mirbooking-sequence.h',
    'mirbooking-target.h',
    'mirbooking-target-site.h']

install_headers(mirbooking_headers, subdir: 'mirbooking-1.0')

mirbooking_sources = [
    'mirbooking-broker.c',
    'mirbooking-error.c',
    'mirbooking-mirna.c',
    'mirbooking-precomputed-score-table.c',
    'mirbooking-score-table.c',
    'mirbooking-sequence.c',
    'mirbooking-target.c']

mirbooking_lib = library('mirbooking-1.0', mirbooking_sources,
                         dependencies: [glib_dep, gobject_dep, gio_dep, m_dep],
                         install: true)

mirbooking_dep = declare_dependency(include_directories: include_directories('.'),
                                    link_with: mirbooking_lib)


if get_option('with_introspection')
    gnome = import('gnome')
    mirbooking_gir = gnome.generate_gir(mirbooking_lib,
                                        sources: mirbooking_headers + mirbooking_sources,
                                        namespace: 'Mirbooking',
                                        nsversion: '1.0',
                                        includes: ['GLib-2.0', 'GObject-2.0', 'Gio-2.0'],
                                        install: true)

    # FIXME: install in 'platlib' see https://github.com/mesonbuild/meson/pull/2196
    install_data('gi/overrides/Mirbooking.py',
                 install_dir: join_paths(import('python3').sysconfig_path('purelib'), 'gi/overrides'))

    if get_option('with_vapi')
        gnome.generate_vapi('mirbooking-1.0',
                            sources: mirbooking_gir[0],
                            packages: ['glib-2.0', 'gobject-2.0', 'gio-2.0'],
                            install: true)
    endif
endif

pkgconfig = import('pkgconfig')
pkgconfig.generate(name: meson.project_name(),
                   description: 'miRBooking simulates the stoichiometric mode of action of microRNAs',
                   version: meson.project_version(),
                   requires: ['glib-2.0', 'gobject-2.0', 'gio-2.0'],
                   libraries: mirbooking_lib,
                   libraries_private: '-lm',
                   filebase: 'mirbooking-1.0',
                   subdirs: 'mirbooking-1.0')
