# disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

add_subdirectory(unit)
add_subdirectory(regression)
add_subdirectory(integration)
if(ENABLE_JSON)
add_subdirectory(examples)
endif()
add_subdirectory(tutorial)
add_subdirectory(surface_rxn)
