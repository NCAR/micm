################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME openmp SOURCES test_openmp.cpp)
if(MICM_ENABLE_LLVM)
  create_standard_test(NAME openmp_jit_solver SOURCES test_openmp_jit_solver.cpp)
endif()