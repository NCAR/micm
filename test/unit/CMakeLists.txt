add_subdirectory(configure)
add_subdirectory(process)
add_subdirectory(solver)
add_subdirectory(system)
add_subdirectory(version)

################################################################################
# Copy test data

add_custom_target(copy_unit_test_configs ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/unit_configs ${CMAKE_BINARY_DIR}/unit_configs)