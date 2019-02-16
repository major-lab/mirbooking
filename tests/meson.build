mirbooking_tests = [
    'score-table',
    'sequence',
    'mirbooking']
foreach name : mirbooking_tests
    test(name, executable(name + '-test', name + '-test.c',
                          dependencies: [glib_dep, gobject_dep, mirbooking_dep, m_dep]))
endforeach

test('integration', mirbooking_bin,
     args: ['--targets', join_paths(meson.current_source_dir(), 'data/gencode.v29.transcripts.fa'),
            '--mirnas', join_paths(meson.current_source_dir(), 'data/mature.fa'),
            '--seed-scores', join_paths(meson.current_source_dir(), '../data/scores-7mer-3mismatch-ending'),
            '--cutoff', '0',
            '--input', join_paths(meson.current_source_dir(), 'data/quantifications.tsv'),
            '--output', '/dev/null'])


if get_option('with_introspection')
    test('gi', python3_bin,
         args: [join_paths(meson.current_source_dir(), 'gi-test.py')],
         env: ['GI_TYPELIB_PATH=@0@'.format(join_paths(meson.current_build_dir(), '../src')),
               'LD_LIBRARY_PATH=@0@'.format(join_paths(meson.current_build_dir(), '../src')),
               'PYTHONPATH=@0@'.format(join_paths(meson.current_source_dir(), '../src'))])
endif
