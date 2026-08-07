################################################################################
# Chapman micro-benchmark
#
# Standalone benchmark that runs the Chapman mechanism at large grid-cell counts
# so per-Solve() work dominates over per-call overhead. Used to check for
# performance regressions in the linear solver and related hot paths.
#
# Not registered as a ctest. Run directly (see scripts/bench_chapman.sh) or
# under callgrind (see scripts/profile_chapman.sh).

add_executable(chapman_bench chapman_bench.cpp)
target_link_libraries(chapman_bench PRIVATE musica::micm)
# Keep full optimization but also emit debug info so callgrind_annotate can
# attribute instruction counts to specific source lines. -g at -O3 has no
# runtime cost.
target_compile_options(chapman_bench PRIVATE $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang,NVHPC>:-g>)
set_target_properties(chapman_bench PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

if(NOT ${MICM_GPU_TYPE} STREQUAL "None")
  target_link_libraries(chapman_bench PRIVATE micm_cuda)
endif()
