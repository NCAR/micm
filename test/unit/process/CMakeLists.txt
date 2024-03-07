################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME process SOURCES test_process.cpp)
create_standard_test(NAME process_set SOURCES test_process_set.cpp)
create_standard_test(NAME arrhenius_rate_constant SOURCES test_arrhenius_rate_constant.cpp)
create_standard_test(NAME branched_rate_constant SOURCES test_branched_rate_constant.cpp)
create_standard_test(NAME surface_rate_constant SOURCES test_surface_rate_constant.cpp)
create_standard_test(NAME ternary_chemical_activation_rate_constant SOURCES test_ternary_chemical_activation_rate_constant.cpp)
create_standard_test(NAME troe_rate_constant SOURCES test_troe_rate_constant.cpp)
create_standard_test(NAME tunneling_rate_constant SOURCES test_tunneling_rate_constant.cpp)
create_standard_test(NAME user_defined_rate_constant SOURCES test_user_defined_rate_constant.cpp)

# GPU tests
if(MICM_ENABLE_CUDA)
  create_standard_test(NAME cuda_process_set SOURCES test_cuda_process_set.cpp LIBRARIES musica::micm_cuda)
endif()

if(MICM_ENABLE_OPENACC)
  create_standard_test(NAME openacc_process_set SOURCES test_openacc_process_set.cpp LIBRARIES musica::micm_openacc)
endif()

if(MICM_ENABLE_LLVM)
  create_standard_test(NAME jit_process_set SOURCES test_jit_process_set.cpp)
  create_standard_test(NAME jit_forcing_calculation SOURCES test_jit_forcing_calculation.cpp forcing_calculation.cpp)
endif()