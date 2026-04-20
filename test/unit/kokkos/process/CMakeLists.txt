add_executable(test_kokkos_process_set test_kokkos_process_set.cpp)
target_link_libraries(test_kokkos_process_set PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_process_set COMMAND test_kokkos_process_set WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
