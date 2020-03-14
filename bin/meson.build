conf = configuration_data()

conf.set('MIRBOOKING_VERSION', meson.project_version())
conf.set('MIRBOOKING_DEFAULT_SEED_SCORES', join_paths(get_option('prefix'), get_option('datadir'), 'mirbooking-@0@'.format(api_version), 'scores-7mer-3mismatch-ending'))
conf.set('MIRBOOKING_DEFAULT_WEE_ET_AL_2012_SUPPLEMENTARY_SCORES', join_paths(get_option('prefix'), get_option('datadir'), 'mirbooking-@0@'.format(api_version), 'scores-4mer'))
conf.set('MIRBOOKING_DEFAULT_YAN_ET_AL_2018_SUPPLEMENTARY_SCORES', join_paths(get_option('prefix'), get_option('datadir'), 'mirbooking-@0@'.format(api_version), 'scores-3mer'))

configure_file(input: 'config.h.in',
    output: 'config.h',
    configuration: conf)
mirbooking_bin = executable('mirbooking', 'mirbooking.c',
                            dependencies: [glib_dep, gobject_dep, gio_dep, gio_unix_dep, mirbooking_dep, m_dep, mpi_dep],
                            install: true)
