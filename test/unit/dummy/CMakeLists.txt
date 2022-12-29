################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME dummy SOURCES test_dummy.cpp)

################################################################################
# OpenMP Tests
if(ENABLE_OPENMP)
  create_standard_test(NAME dummy_openmp SOURCES test_dummy_openmp.cpp)
endif()

################################################################################
# MPI Tests
if(ENABLE_MPI)
  create_standard_test(NAME dummy_mpi SOURCES test_dummy_mpi.cpp)
endif()