################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

create_standard_test(NAME phase SOURCES test_phase.cpp)
create_standard_test(NAME species SOURCES test_species.cpp)
create_standard_test(NAME system SOURCES test_system.cpp)
