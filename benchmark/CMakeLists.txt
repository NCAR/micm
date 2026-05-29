# disable clang tidy for benchmark files
set(CMAKE_CXX_CLANG_TIDY "")

################################################################################
# DAE / ODE benchmarks (built only when MICM_ENABLE_BENCHMARKS=ON)

add_executable(robertson_dae robertson_dae.cpp)
target_link_libraries(robertson_dae PUBLIC musica::micm)
