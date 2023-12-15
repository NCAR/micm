# disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

add_executable(micmDriver example.cpp)
target_link_libraries(micmDriver PUBLIC musica::micm )

set_target_properties(micmDriver PROPERTIES OUTPUT_NAME "micm")

################################################################################
# Run each example configuration as a test

add_test(NAME carbon_bond_5_example
         COMMAND micm ${CMAKE_BINARY_DIR}/examples/configs/carbon_bond_5 ${CMAKE_BINARY_DIR}/examples/configs/carbon_bond_5/initial_conditions.csv
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME chapman_example
         COMMAND micm ${CMAKE_BINARY_DIR}/examples/configs/chapman ${CMAKE_BINARY_DIR}/examples/configs/chapman/initial_conditions.csv
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME robertson_example
         COMMAND micm ${CMAKE_BINARY_DIR}/examples/configs/robertson ${CMAKE_BINARY_DIR}/examples/configs/robertson/initial_conditions.csv
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME ts1_example
         COMMAND micm ${CMAKE_BINARY_DIR}/examples/configs/TS1 ${CMAKE_BINARY_DIR}/examples/configs/TS1/initial_conditions.csv
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})