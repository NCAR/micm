################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME chapman_ode_solver SOURCES test_chapman_ode_solver.cpp)
create_standard_test(NAME chapman_ode_solver_builder SOURCES test_chapman_ode_solver_builder.cpp)
create_standard_test(NAME rosenbrock SOURCES test_rosenbrock.cpp)