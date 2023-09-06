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

# GPU tests
if(ENABLE_CUDA)
  create_standard_test(NAME cuda_lu_decomposition SOURCES test_cuda_lu_decomposition.cpp LIBRARIES musica::micm_cuda)
endif()

if(ENABLE_OPENACC)
  create_standard_test(NAME openacc_process_set SOURCES test_openacc_process_set.cpp LIBRARIES musica::micm_openacc)
endif()

if(ENABLE_LLVM)
  create_standard_test(NAME jit_process_set SOURCES test_jit_process_set.cpp)
endif()