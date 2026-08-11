################################################################################
# Tests
add_executable(test_kokkos_constraint_set test_kokkos_constraint_set.cpp)
target_link_libraries(test_kokkos_constraint_set PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_constraint_set COMMAND test_kokkos_constraint_set WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(type)