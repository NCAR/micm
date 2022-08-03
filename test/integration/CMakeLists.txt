################################################################################
# micm integration tests

add_executable(test_c_host c_host.c)
target_link_libraries(test_c_host PUBLIC micmlib)
add_test(NAME c_host COMMAND test_c_host)

################################################################################
