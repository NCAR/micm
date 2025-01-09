################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME chapman_ode_solver SOURCES test_chapman_ode_solver.cpp)
create_standard_test(NAME linear_solver SOURCES test_linear_solver.cpp)
create_standard_test(NAME linear_solver_in_place SOURCES test_linear_solver_in_place.cpp)
create_standard_test(NAME lu_decomposition_doolittle SOURCES test_lu_decomposition_doolittle.cpp)
create_standard_test(NAME lu_decomposition_mozart SOURCES test_lu_decomposition_mozart.cpp)
create_standard_test(NAME lu_decomposition_doolittle_in_place SOURCES test_lu_decomposition_doolittle_in_place.cpp)
create_standard_test(NAME lu_decomposition_mozart_in_place SOURCES test_lu_decomposition_mozart_in_place.cpp)
create_standard_test(NAME rosenbrock SOURCES test_rosenbrock.cpp)
create_standard_test(NAME state SOURCES test_state.cpp)
create_standard_test(NAME backward_euler SOURCES test_backward_euler.cpp)
create_standard_test(NAME solver_builder SOURCES test_solver_builder.cpp)