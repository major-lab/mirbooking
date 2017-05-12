install_headers(['mirbooking.h',
                 'mirbooking-mirna.h',
                 'mirbooking-score-table.h',
                 'mirbooking-sequence.h',
                 'mirbooking-target.h'], subdir: 'mirbooking-1.0')

mirbooking_sources = [
    'mirbooking.c',
    'mirbooking-mirna.c',
    'mirbooking-score-table.c',
    'mirbooking-sequence.c',
    'mirbooking-target.c']

mirbooking_lib = library('mirbooking', mirbooking_sources,
                         dependencies: [glib_dep, gobject_dep],
                         install: true)

mirbooking_dep = declare_dependency(include_directories: include_directories('.'),
                                    link_with: mirbooking_lib)
