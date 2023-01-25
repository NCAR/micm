################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME cpp_fortran_interface SOURCES interface.cpp)
create_standard_test(NAME fortran_can_call_cpp SOURCES interface.F90)
