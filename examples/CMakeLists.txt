# disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

add_executable(micmDriver example.cpp)
target_link_libraries(micmDriver PUBLIC musica::micm )

set_target_properties(micmDriver PROPERTIES OUTPUT_NAME "micm")