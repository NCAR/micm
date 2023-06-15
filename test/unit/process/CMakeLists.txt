################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME arrhenius_rate_constant SOURCES test_arrhenius_rate_constant.cpp)
create_standard_test(NAME photolysis_rate_constant SOURCES test_photolysis_rate_constant.cpp)
create_standard_test(NAME troe_rate_constant SOURCES test_troe_rate_constant.cpp)
create_standard_test(NAME process_set SOURCES test_process_set.cpp)

# GPU tests
if(ENABLE_CUDA)
  create_standard_test(NAME cuda_process_set SOURCES test_cuda_process_set.cpp)
endif()

if(ENABLE_OPENACC)
  create_standard_test(NAME openacc_process_set SOURCES test_openacc_process_set.cpp)
endif()