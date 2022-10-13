################################################################################
# Link micmlib to a test and add it to the suite

macro(add_std_test test_name)
  target_include_directories(${test_name} PUBLIC ${CMAKE_BINARY_DIR}/src)
  target_link_libraries(${test_name} micmlib)
  add_test(NAME ${test_name} COMMAND ${test_name})
endmacro(add_std_test)

################################################################################
# MICM performance tests

add_executable(test_performance driver.F90)
add_std_test(test_performance)

################################################################################

