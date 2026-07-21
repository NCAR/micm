# disable clang tidy for benchmark files
set(CMAKE_CXX_CLANG_TIDY "")

################################################################################
# DAE / ODE benchmarks (built only when MICM_ENABLE_BENCHMARKS=ON)

add_executable(robertson_dae robertson_dae.cpp)
target_link_libraries(robertson_dae PUBLIC musica::micm)
target_compile_definitions(robertson_dae PRIVATE MICM_BENCHMARK_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

add_executable(equilibrium_efficiency equilibrium_efficiency.cpp)
target_link_libraries(equilibrium_efficiency PUBLIC musica::micm)

add_executable(tropospheric_dae tropospheric_dae.cpp)
target_link_libraries(tropospheric_dae PUBLIC musica::micm)
target_compile_definitions(tropospheric_dae PRIVATE MICM_BENCHMARK_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")

add_executable(dae_convergence_experiments dae_convergence_experiments.cpp)
target_link_libraries(dae_convergence_experiments PUBLIC musica::micm)

add_executable(test_robertson_qssa test_robertson_qssa.cpp)
target_link_libraries(test_robertson_qssa PUBLIC musica::micm GTest::gtest_main)
add_test(NAME test_robertson_qssa COMMAND test_robertson_qssa)

add_executable(robertson_conservation_dae robertson_conservation_dae.cpp)
target_link_libraries(robertson_conservation_dae PUBLIC musica::micm)
target_compile_definitions(robertson_conservation_dae PRIVATE MICM_BENCHMARK_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
