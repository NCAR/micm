################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

add_executable(regression_test_rosebrock_p_force regression_test_p_force.cpp)
target_link_libraries(regression_test_rosebrock_p_force musica::micm micmlib GTest::gtest_main)
add_test(NAME rosenbrock_p_force
          COMMAND regression_test_rosebrock_p_force ""
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
