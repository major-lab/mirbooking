install_headers(['mirbooking.h',
                 'mirbooking-mirna.h',
                 'mirbooking-score-table.h',
                 'mirbooking-sequence.h',
                 'mirbooking-target.h'], subdir: 'mirbooking-1.0')

mirbooking_sources = [
    'mirbooking-broker.c',
    'mirbooking-mirna.c',
    'mirbooking-score-table.c',
    'mirbooking-sequence.c',
    'mirbooking-target.c']

mirbooking_lib = library('mirbooking-1.0', mirbooking_sources,
                         dependencies: [glib_dep, gobject_dep, m_dep],
                         install: true)

mirbooking_dep = declare_dependency(include_directories: include_directories('.'),
                                    link_with: mirbooking_lib)

pkgconfig = import('pkgconfig')
pkgconfig.generate(name: meson.project_name(),
                   description: 'miRBooking simulates the stoichiometric mode of action of microRNAs',
                   version: meson.project_version(),
                   requires: ['glib-2.0'],
                   libraries: mirbooking_lib,
                   libraries_private: '-lm',
                   filebase: 'mirbooking-1.0',
                   subdirs: 'mirbooking-1.0')
