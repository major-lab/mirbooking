executable('mirbooking', 'mirbooking.c',
           dependencies: [glib_dep, gobject_dep, gio_dep, mirbooking_dep, m_dep],
           install: true)

install_data(['mirbooking-aggregate',
              'mirbooking-calibrate',
              'mirbooking-extract-cds-regions'], install_dir: 'bin')
