// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Standalone benchmark for measuring solver performance.
//
// Runs one mechanism at a large grid-cell count, so per-Solve() work dominates
// over per-call overhead.
//
// Usage: micm_bench [num_cells] [num_steps] [dt_seconds] [backend] [matrix] [lu_type] [lu_algorithm] [mechanism]
//   backend:      "cpu" (default) or "gpu"
//   matrix:       "standard" (default) or "vector1|2|4|8|128"
//   lu_type:      "in-place" (default) or "separate"
//   lu_algorithm: "mozart" (default) or "doolittle"
//   mechanism:    "chapman" (default) or "ts1"
//
//   Note that "gpu" backend must use "mozart"/"in-place" LU and "vector1|2|4|8|128" matrix
//
// Pass "list" as the first argument to print every valid combination.
//
// Prints a single line with configuration and elapsed wall time (ms).

#include "bench_harness.hpp"
#include "chapman_mechanism.hpp"
#include "ts1_mechanism.hpp"

#ifdef MICM_USE_CUDA
  #include <micm/cuda/util/cuda_util.cuh>
#endif

#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

namespace
{
  bench::Registry BuildRegistry()
  {
    bench::Registry registry;
    bench::Register<bench::Chapman>(registry, bench::VectorWidths{});
    bench::Register<bench::Ts1>(registry, bench::VectorWidths{});
    return registry;
  }

  constexpr const char* kUsage =
      " [num_cells] [num_steps] [dt_seconds] [backend] [matrix] [lu_type] [lu_algorithm] [mechanism]";
}  // namespace

int main(int argc, char** argv)
{
  const auto registry = BuildRegistry();

  if (argc > 1 && std::string(argv[1]) == "list")
  {
    for (const auto& [key, runner] : registry)
    {
      std::cout << key << std::endl;
    }
    return 0;
  }

  // std::stoul and std::stod throw on a non-numeric argument, so parse inside a
  // try block and report the usage rather than call std::terminate.
  bench::Config config{};
  try
  {
    config = bench::Config{
      .num_cells_ = (argc > 1) ? static_cast<micm::Index>(std::stoul(argv[1])) : micm::Index{ 10000 },
      .num_steps_ = (argc > 2) ? static_cast<micm::Index>(std::stoul(argv[2])) : micm::Index{ 100 },
      .dt_ = (argc > 3) ? static_cast<micm::Real>(std::stod(argv[3])) : static_cast<micm::Real>(30.0),
      .backend_ = (argc > 4) ? argv[4] : "cpu",
      .matrix_ = (argc > 5) ? argv[5] : "standard",
      .lu_matrix_ = (argc > 6) ? argv[6] : "in-place",
      .lu_ = (argc > 7) ? argv[7] : "mozart",
      .mechanism_ = (argc > 8) ? argv[8] : "chapman",
    };
  }
  catch (const std::exception& e)
  {
    std::cout << "Could not read the arguments: " << e.what() << std::endl;
    std::cout << "Usage: " << argv[0] << kUsage << std::endl;
    return 1;
  }

  const auto entry = registry.find(bench::Key(config));
  if (entry == registry.end())
  {
    std::cout << "Invalid option combination: " << bench::Key(config) << std::endl;
    std::cout << "Run '" << argv[0] << " list' to print every valid combination." << std::endl;
    return 1;
  }

  double elapsed_ms = -1.0;
  try
  {
    elapsed_ms = entry->second(config);
  }
  catch (const std::exception& e)
  {
    std::cout << "Benchmark failed: " << e.what() << std::endl;
    return 1;
  }

#ifdef MICM_USE_CUDA
  micm::cuda::CudaStreamSingleton::GetInstance().CleanUp();
#endif

  // Fixed notation, always. The default format switches to scientific above
  // 1e6 ms, and the "elapsed_ms=([0-9.]+)" pattern in bench_micm.sh would then
  // read 1.23457e+06 as 1.23457.
  std::cout << "mechanism=" << config.mechanism_ << " backend=" << config.backend_ << " matrix_type=" << config.matrix_
            << " lu=" << config.lu_ << "/" << config.lu_matrix_ << " cells=" << config.num_cells_
            << " steps=" << config.num_steps_ << " dt=" << config.dt_ << std::fixed << std::setprecision(3)
            << " elapsed_ms=" << elapsed_ms << " ms_per_step=" << elapsed_ms / static_cast<double>(config.num_steps_)
            << std::endl;
  return 0;
}
