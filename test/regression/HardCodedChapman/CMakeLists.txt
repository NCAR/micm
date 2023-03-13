################################################################################
# Test utilities

include(test_util)

################################################################################
# Tests

add_executable(regression_test_arrhenius regression_test_arrhenius.cpp rate_constants.F90)
target_link_libraries(regression_test_arrhenius musica::micm micmlib GTest::gtest_main)
add_test(NAME arrhenius_regression
          COMMAND regression_test_arrhenius ""
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(regression_test_name_retrieval regression_test_name_retrieval.cpp name_retrieval.F90)
target_link_libraries(regression_test_name_retrieval musica::micm micmlib GTest::gtest_main)
add_test(NAME name_regression
          COMMAND regression_test_name_retrieval ""
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(regression_test_p_force regression_test_p_force.cpp p_force.F90)
target_link_libraries(regression_test_p_force musica::micm micmlib GTest::gtest_main)
add_test(NAME p_force
          COMMAND regression_test_p_force ""
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(regression_test_factored_alpha_minus_jac regression_test_factored_alpha_minus_jac.cpp factored_alpha_minus_jac.F90)
target_link_libraries(regression_test_factored_alpha_minus_jac musica::micm micmlib GTest::gtest_main)
add_test(NAME factored_alpha_minus_jac
          COMMAND regression_test_factored_alpha_minus_jac ""
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(regression_test_dforce_dy_times_vector regresstion_test_dforce_dy_times_vector.cpp dforce_dy_times_vector.F90)
target_link_libraries(regression_test_dforce_dy_times_vector musica::micm micmlib GTest::gtest_main)
add_test(NAME dforce_dy_times_vector
          COMMAND regression_test_dforce_dy_times_vector ""
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(regression_test_kinetic_solve regression_test_kinetic_solve.cpp kinetic_solve.F90)
target_link_libraries(regression_test_kinetic_solve musica::micm micmlib GTest::gtest_main)
add_test(NAME kinetic_solve
          COMMAND regression_test_kinetic_solve ""
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(regression_test_dforce_dy regression_test_dforce_dy.cpp dforce_dy.F90)
target_link_libraries(regression_test_dforce_dy musica::micm micmlib GTest::gtest_main)
add_test(NAME dforce_dy
          COMMAND regression_test_dforce_dy ""
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(regression_test_solve regression_test_solve.cpp solve.F90)
target_link_libraries(regression_test_solve musica::micm micmlib GTest::gtest_main)
add_test(NAME solve
          COMMAND regression_test_solve ""
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
