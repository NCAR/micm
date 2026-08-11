################################################################################
# Tests
add_executable(test_kokkos_equilibrium test_kokkos_equilibrium.cpp)
target_link_libraries(test_kokkos_equilibrium PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_equilibrium COMMAND test_kokkos_equilibrium WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_kokkos_linear test_kokkos_linear.cpp)
target_link_libraries(test_kokkos_linear PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_linear COMMAND test_kokkos_linear WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
