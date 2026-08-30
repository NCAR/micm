# Clang-tidy is disabled for benchmark programs: they are measurement harnesses,
# not library code, and the repository's tidy profile is tuned for the latter.
set(CMAKE_CXX_CLANG_TIDY "")

################################################################################
# DAE constraint-initialization benchmarks (built only when
# MICM_ENABLE_DAE_BENCHMARKS=ON). These register no ctest: they emit CSV for the
# evidence in docs/superpowers/notes/, and are run by hand.

add_executable(dae_init_cold_start dae_init_cold_start.cpp)
target_link_libraries(dae_init_cold_start PUBLIC musica::micm)
