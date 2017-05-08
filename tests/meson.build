mirbooking_tests = [
    'score-table',
    'sequence',
    'target']
foreach name : mirbooking_tests
    test(name, executable(name + '-test', name + '-test.c',
                          dependencies: [glib_dep, gobject_dep, mirbooking_dep]))
endforeach
