mirbooking_headers = [
    'mirbooking.h',
    'mirbooking-error.h',
    'mirbooking-broker.h',
    'mirbooking-broker-sparse-solver.h',
    'mirbooking-broker-output-format.h',
    'mirbooking-mirna.h',
    'mirbooking-occupant.h',
    'mirbooking-default-score-table.h',
    'mirbooking-score-table.h',
    'mirbooking-sequence.h',
    'mirbooking-target.h',
    'mirbooking-target-site.h',
    'mirbooking-utils.h']

install_headers(mirbooking_headers, subdir: 'mirbooking-@0@'.format(api_version))

mirbooking_sources = [
    'mirbooking-broker.c',
    'mirbooking-broker-sparse-solver.c',
    'mirbooking-broker-output-format.c',
    'mirbooking-error.c',
    'mirbooking-mirna.c',
    'mirbooking-occupant.c',
    'mirbooking-default-score-table.c',
    'mirbooking-score-table.c',
    'mirbooking-sequence.c',
    'mirbooking-target.c',
    'mirbooking-utils.c']

mirbooking_lib = library('mirbooking-@0@'.format(api_version), mirbooking_sources,
                         dependencies: [glib_dep, gobject_dep, gio_dep, m_dep, pb_dep, openmp_dep, mpi_dep, odeint_dep, sparse_dep, cblas_dep],
                         install: true)

mirbooking_dep = declare_dependency(include_directories: include_directories('.'),
                                    link_with: mirbooking_lib)


if get_option('with_introspection')
    gnome = import('gnome')
    mirbooking_gir = gnome.generate_gir(mirbooking_lib,
                                        sources: mirbooking_headers + mirbooking_sources,
                                        namespace: 'Mirbooking',
                                        nsversion: api_version,
                                        includes: ['GLib-2.0', 'GObject-2.0', 'Gio-2.0'],
                                        install: true)

    if python3.found()
        python3.install_sources('gi/overrides/Mirbooking.py',
                                subdir: 'gi/overrides',
                                pure: false)
    endif

    if get_option('with_vapi')
        gnome.generate_vapi('mirbooking-2.0',
                            sources: mirbooking_gir[0],
                            packages: ['glib-2.0', 'gobject-2.0', 'gio-2.0'],
                            install: true)
    endif
endif

pkgconfig = import('pkgconfig')
pkgconfig.generate(mirbooking_lib,
                   name: meson.project_name(),
                   description: 'miRBooking simulates the stoichiometric mode of action of microRNAs',
                   version: meson.project_version(),
                   filebase: 'mirbooking-@0@'.format(api_version),
                   subdirs: 'mirbooking-@0@'.format(api_version))
