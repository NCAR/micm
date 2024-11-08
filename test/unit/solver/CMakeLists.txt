################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME chapman_ode_solver SOURCES test_chapman_ode_solver.cpp)
create_standard_test(NAME linear_solver SOURCES test_linear_solver.cpp)
create_standard_test(NAME doolittle_lu_decomposition SOURCES test_doolittle_lu_decomposition.cpp)
create_standard_test(NAME rosenbrock SOURCES test_rosenbrock.cpp)
create_standard_test(NAME state SOURCES test_state.cpp)
create_standard_test(NAME backward_euler SOURCES test_backward_euler.cpp)
create_standard_test(NAME solver_builder SOURCES test_solver_builder.cpp)