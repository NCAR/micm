################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

if(ENABLE_JSON)
  create_standard_test(NAME surface_rxn SOURCES test_surface_rxn.cpp)
endif()

################################################################################

################################################################################
# Copy test data

add_custom_target(copy_kpp_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/configs ${CMAKE_BINARY_DIR}/configs)

################################################################################
