# disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

add_executable(micmprofile profile.cpp initial_conditions.hpp)
target_link_libraries(micmprofile PUBLIC musica::micm )

add_test(NAME ts1_profiling
         COMMAND micmprofile ${CMAKE_BINARY_DIR}//TS1 ${CMAKE_BINARY_DIR}/examples/configs/TS1/initial_conditions.csv
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})