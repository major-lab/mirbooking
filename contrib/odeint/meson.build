odeint_lib = static_library('odeint', 'odeint.c', dependencies: [m_dep, cblas_dep])
odeint_dep = declare_dependency(include_directories: include_directories('.'),
                                link_with: odeint_lib)
test('odeint', executable('odeint-test', 'odeint-test.c', dependencies: odeint_dep))
