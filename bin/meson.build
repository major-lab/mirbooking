mirbooking_bin = executable('mirbooking', 'mirbooking.c',
                            dependencies: [glib_dep, gobject_dep, gio_dep, mirbooking_dep, m_dep],
                            install: true)

install_data(['mirbooking-calibrate',
              'mirbooking-extract-cds-regions',
              'mirbooking-generate-score-table'], install_dir: 'bin')
