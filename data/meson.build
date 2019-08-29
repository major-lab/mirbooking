install_data(['scores-7mer',
              'scores-7mer-1mismatch-ending',
              'scores-7mer-2mismatch-ending',
              'scores-7mer-3mismatch-ending',
              'scores-4mer',
              'scores-3mer'],
             install_dir: 'share/mirbooking-@0@'.format(api_version))

help2man_bin = find_program('help2man', required: false)
if help2man_bin.found()
    custom_target('mirbooking man page', command: [help2man_bin, '--name', 'simulate the stochiometric mode of action of miRNAs', '--version-string', meson.project_version(), '--no-info', '--output', '@OUTPUT@', '@INPUT@'],
                  input: mirbooking_bin,
                  output: 'mirbooking.1',
                  install: true,
                  install_dir: 'share/man/man1')
endif
