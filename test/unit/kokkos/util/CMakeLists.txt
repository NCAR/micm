add_executable(test_kokkos_dense_matrix test_kokkos_dense_matrix.cpp)
target_link_libraries(test_kokkos_dense_matrix PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_dense_matrix COMMAND test_kokkos_dense_matrix WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(test_kokkos_sparse_matrix test_kokkos_sparse_matrix.cpp)
target_link_libraries(test_kokkos_sparse_matrix PUBLIC musica::micm musica::micm_kokkos GTest::gtest)
add_test(NAME kokkos_sparse_matrix COMMAND test_kokkos_sparse_matrix WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
