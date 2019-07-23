mirbooking_tests = [
    'score-table',
    'sequence',
    'mirbooking']
foreach name : mirbooking_tests
    test(name, executable(name + '-test', name + '-test.c',
                          dependencies: [glib_dep, gobject_dep, mirbooking_dep, m_dep]),
                          env: ['G_TEST_SRCDIR=@0@'.format( meson.current_source_dir())])
endforeach

mirbooking_integration_tests = [
    ['superlu',     ['--sparse-solver', 'superlu'],     get_option('with_superlu')],
    ['umfpack',     ['--sparse-solver', 'umfpack'],     get_option('with_umfpack')],
    ['mkl-dss',     ['--sparse-solver', 'mkl-dss'],     get_option('with_mkl_dss')],
    ['mkl-cluster', ['--sparse-solver', 'mkl-cluster'], get_option('with_mkl_cluster')],
    ['cusolver',    ['--sparse-solver', 'cusolver'],    get_option('with_cusolver')],
    ['pardiso',     ['--sparse-solver', 'pardiso'],     get_option('with_pardiso')],
    ['refseq',       ['--targets', join_paths(meson.current_source_dir(), 'data/GCF_000001405.39_GRCh38.p13_rna.fna'), '--input', join_paths(meson.current_source_dir(), 'data/quantifications-refseq.tsv')], true]
]

mirbooking_integration_test_args = [
    '--targets', join_paths(meson.current_source_dir(), 'data/gencode.v29.transcripts.fa'),
    '--mirnas', join_paths(meson.current_source_dir(), 'data/mature.fa'),
    '--seed-scores', join_paths(meson.current_source_dir(), '../data/scores-7mer-3mismatch-ending'),
    '--supplementary-scores', join_paths(meson.current_source_dir(), '../data/scores-3mer'),
    '--cutoff', '0',
    '--input', join_paths(meson.current_source_dir(), 'data/quantifications.tsv'),
    '--output', '/dev/null']

foreach t : mirbooking_integration_tests
    if t[2]
        test('integration-@0@'.format(t[0]), mirbooking_bin, args: mirbooking_integration_test_args + t[1])
    endif
endforeach

if get_option('with_introspection')
    test('gi', python3,
         args: [join_paths(meson.current_source_dir(), 'gi-test.py')],
         env: ['GI_TYPELIB_PATH=@0@'.format(join_paths(meson.current_build_dir(), '../src')),
               'LD_LIBRARY_PATH=@0@'.format(join_paths(meson.current_build_dir(), '../src')),
               'PYTHONPATH=@0@'.format(join_paths(meson.current_source_dir(), '../src'))])
endif
