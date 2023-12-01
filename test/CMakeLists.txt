# disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

add_subdirectory(unit)
add_subdirectory(regression)
add_subdirectory(integration)
add_subdirectory(tutorial)
