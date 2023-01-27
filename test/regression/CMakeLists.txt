################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

# compile the old micm library
add_subdirectory(old_micm)

add_executable(regression_test_arrhenius test_arrhenius_regression.cpp rate_constants.F90)
target_link_libraries(regression_test_arrhenius musica::micm micmlib GTest::gtest_main)

add_test(NAME arrhenius_regression
          COMMAND regression_test_arrhenius ""
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})