################################################################################
# DAE / ODE benchmarks (built only when MICM_ENABLE_BENCHMARKS=ON)

add_executable(robertson_dae robertson_dae.cpp)
target_link_libraries(robertson_dae PUBLIC musica::micm)
target_compile_features(robertson_dae PUBLIC cxx_std_20)
