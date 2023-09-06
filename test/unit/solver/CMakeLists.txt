################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME chapman_ode_solver SOURCES test_chapman_ode_solver.cpp)
create_standard_test(NAME linear_solver SOURCES test_linear_solver.cpp)
create_standard_test(NAME lu_decomposition SOURCES test_lu_decomposition.cpp)
create_standard_test(NAME rosenbrock SOURCES test_rosenbrock.cpp)
create_standard_test(NAME state SOURCES test_state.cpp)

if(ENABLE_LLVM)
  create_standard_test(NAME jit_rosenbrock SOURCES test_jit_rosenbrock.cpp)
endif()