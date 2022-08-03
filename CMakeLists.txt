cmake_minimum_required(VERSION 3.13)
project(micm)
enable_language(C)
enable_language(CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

################################################################################
# micm headers

include_directories(${PROJECT_SOURCE_DIR}/include)

################################################################################
# micm source

add_subdirectory(src)
include_directories(${CMAKE_BINARY_DIR}/src)

################################################################################
# micm tests

enable_testing()
add_subdirectory(test)

################################################################################
