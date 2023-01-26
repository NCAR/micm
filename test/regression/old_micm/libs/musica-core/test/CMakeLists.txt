################################################################################
# Link musicalib to a test and add it to the suite

macro(add_std_test test_name)
  target_include_directories(${test_name} PUBLIC ${CMAKE_BINARY_DIR}/src)
  target_link_libraries(${test_name} musicalib)
  add_test(NAME ${test_name} COMMAND ${test_name})
endmacro(add_std_test)

################################################################################
# Link musicalib to a test and add it to the suite as a bash script

macro(add_std_test_script test_name script_path)
  target_include_directories(${test_name} PUBLIC ${CMAKE_BINARY_DIR}/src)
  target_link_libraries(${test_name} musicalib)
  add_test(NAME ${test_name} COMMAND ${script_path})
endmacro(add_std_test_script)

################################################################################
# Copy test data

add_custom_target(copy_unit_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/test)

################################################################################
# MusicBox tests

add_executable(config config.F90)
add_std_test(config)

add_executable(domain_cell domains/cell.F90)
add_std_test(domain_cell)

add_executable(grid grid.F90)
add_std_test(grid)

add_executable(initial_conditions initial_conditions.F90
                                  mock_domain.F90
                                  mock_domain_iterator.F90
                                  mock_domain_state.F90
                                  mock_domain_state_accessor.F90
                                  mock_domain_state_mutator.F90)
add_std_test(initial_conditions)

add_executable(interpolator interpolator.F90)
add_std_test(interpolator)

add_executable(interpolator_linear_1D interpolator_strategies/linear_1D.F90)
add_std_test(interpolator_linear_1D)

add_executable(io_netcdf io/netcdf.F90)
add_std_test(io_netcdf)

add_executable(lookup_axis lookup/axis.F90)
add_std_test_script(lookup_axis lookup/axis.sh)

add_executable(lookup_2D_axis lookup/2D_axis.F90)
add_std_test(lookup_2D_axis)

add_executable(util_array util/array.F90)
add_std_test(util_array)

add_executable(util_assert util/assert.F90)
add_std_test_script(util_assert util/assert.sh)

add_executable(util_convert util/convert.F90)
add_std_test(util_convert)

add_executable(util_data_type util/data_type.F90)
add_std_test(util_data_type)

add_executable(util_datetime util/datetime.F90)
add_std_test(util_datetime)

add_executable(util_logger util/logger.F90)
add_std_test(util_logger)

add_executable(util_math util/math.F90)
add_std_test_script(util_math util/math.sh)

add_executable(util_property util/property.F90)
add_std_test(util_property)

add_executable(util_property_set util/property_set.F90)
add_std_test_script(util_property_set util/property_set.sh)

add_executable(util_string util/string.F90)
add_std_test_script(util_string util/string.sh)

################################################################################
