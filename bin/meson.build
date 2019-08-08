mirbooking_bin = executable('mirbooking', 'mirbooking.c',
                            dependencies: [glib_dep, gobject_dep, gio_dep, mirbooking_dep, m_dep, mpi_dep],
                            install: true)

executable('mirbooking-generate-score-table', 'mirbooking-generate-score-table.c',
           dependencies: [glib_dep, gio_dep, m_dep, openmp_dep, sparse_dep],
           install: true)
