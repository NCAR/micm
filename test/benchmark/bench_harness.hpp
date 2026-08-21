// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Shared harness for the standalone solver benchmarks.
//
// The benchmark runs one mechanism against one solver configuration, and the
// two choices are independent:
//
//   * A mechanism is a type with a builder-generic Build() and its own
//     InitState(). See chapman_mechanism.hpp and ts1_mechanism.hpp.
//   * A solver configuration is a matrix ordering, an LU storage choice, an LU
//     algorithm, and a backend.
//
// Register() expands the cross product into a lookup table.
#pragma once

#include <micm/CPU.hpp>
#ifdef MICM_USE_CUDA
  #include <micm/GPU.hpp>
#endif
#ifdef MICM_USE_KOKKOS
  #include <micm/Kokkos.hpp>
#endif

#include <chrono>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

// Callgrind client-request macros. Without valgrind-devel, or outside valgrind,
// these are no-ops.
//
// The profile script passes --collect-atstart=no --instr-atstart=no, and we
// toggle collection on only around the timed Solve() loop, so the reported
// count covers Rosenbrock Solve() and everything it calls.
#if __has_include(<valgrind/callgrind.h>)
  #include <valgrind/callgrind.h>
#else
  #define CALLGRIND_START_INSTRUMENTATION \
    do                                    \
    {                                     \
    } while (0)
  #define CALLGRIND_STOP_INSTRUMENTATION \
    do                                   \
    {                                    \
    } while (0)
  #define CALLGRIND_TOGGLE_COLLECT \
    do                             \
    {                              \
    } while (0)
  #define CALLGRIND_ZERO_STATS \
    do                         \
    {                          \
    } while (0)
#endif

namespace bench
{
  /// @brief One benchmark run, as selected on the command line.
  struct Config
  {
    // Declared in command-line order, so the initializer in main() reads the
    // same way the arguments are passed.
    micm::Index num_cells_;
    micm::Index num_steps_;
    micm::Real dt_;
    std::string backend_;
    std::string matrix_;
    std::string lu_matrix_;
    std::string lu_;
    std::string mechanism_;
  };

  /// @brief Look-up key for one mechanism and one solver configuration.
  inline std::string Key(
      const std::string& mechanism,
      const std::string& backend,
      const std::string& matrix,
      const std::string& lu_matrix,
      const std::string& lu)
  {
    return mechanism + "/" + backend + "/" + matrix + "/" + lu_matrix + "/" + lu;
  }

  inline std::string Key(const Config& config)
  {
    return Key(config.mechanism_, config.backend_, config.matrix_, config.lu_matrix_, config.lu_);
  }

  /// @brief Throw when a solve did not converge.
  ///
  /// A failed solve takes a different number of internal steps, and a NaN in the
  /// state persists, so the measurement would track the failure, not the hot path.
  inline void RequireConverged(micm::SolverState state, const Config& config, const std::string& when)
  {
    if (state == micm::SolverState::Converged)
    {
      return;
    }
    throw std::runtime_error(
        "the " + when + " solve did not converge (" + micm::SolverStateToString(state) + ") for " + Key(config) +
        " at dt=" + std::to_string(config.dt_) + "; reduce dt or correct the initial conditions");
  }

  /// @brief Time the Solve() loop and return the elapsed wall time in ms.
  template<class Solver>
  double TimeSolve(Solver& solver, typename Solver::StatePolicyType& state, const Config& config)
  {
    // Warm-up (first call primes caches / any lazy init).
    auto warmup = solver.Solve(config.dt_, state);
    RequireConverged(warmup.state_, config, "warm-up");

    // Turn callgrind instrumentation on just for the timed loop. When not
    // running under callgrind these macros are no-ops.
    CALLGRIND_ZERO_STATS;
    CALLGRIND_START_INSTRUMENTATION;
    CALLGRIND_TOGGLE_COLLECT;
    auto t0 = std::chrono::steady_clock::now();
    for (micm::Index s = 0; s < config.num_steps_; ++s)
    {
      [[maybe_unused]] auto r = solver.Solve(config.dt_, state);
    }
    auto t1 = std::chrono::steady_clock::now();
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;

    // One more solve, after collection stops, confirms the state stayed healthy.
    // Checking here rather than inside the loop keeps the count measuring Solve().
    auto verify = solver.Solve(config.dt_, state);
    RequireConverged(verify.state_, config, "post-loop");

    return std::chrono::duration<double, std::milli>(t1 - t0).count();
  }

  /// @brief The single benchmark body. Every registered combination runs through here.
  template<class Mechanism, class SolverBuilderPolicy>
  double RunCase(const Config& config)
  {
    auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
    auto solver = Mechanism::Build(SolverBuilderPolicy(options));
    auto state = solver.GetState(config.num_cells_);
    Mechanism::InitState(state, config.num_cells_);
    // For Kokkos, copy inputs to device (no-op for CPU solvers)
    state.variables_.CopyToDevice();
    state.custom_rate_parameters_.CopyToDevice();
    state.conditions_.CopyToDevice();
    // For Kokkos, UpdateStateParameters runs on device with copied over custom_rate_parameters_
    solver.UpdateStateParameters(state);
    // For CUDA, copy the information to the device
    if constexpr (requires { state.SyncInputsToDevice(); })
    {
      state.SyncInputsToDevice();
    }
    return TimeSolve(solver, state, config);
  }

  using Runner = double (*)(const Config&);
  using Registry = std::map<std::string, Runner>;

  using StandardDense = micm::Matrix<micm::Real>;
  template<micm::Index L>
  using VectorDense = micm::VectorMatrix<micm::Real, L>;
  using StandardSparse = micm::SparseMatrix<micm::Real, micm::SparseMatrixStandardOrdering>;
  template<micm::Index L>
  using VectorSparse = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>;

  template<class SM>
  using DooLU = micm::LuDecompositionDoolittle<SM>;
  template<class SM>
  using DooLUInPlace = micm::LuDecompositionDoolittleInPlace<SM>;
  template<class SM>
  using MozLU = micm::LuDecompositionMozart<SM>;
  template<class SM>
  using MozLUInPlace = micm::LuDecompositionMozartInPlace<SM>;

  template<class DM, class SM, class LU>
  using CpuRosen = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters, DM, SM, LU>;
  template<class DM, class SM, class LU>
  using CpuRosenInPlace = micm::CpuSolverBuilderInPlace<micm::RosenbrockSolverParameters, DM, SM, LU>;
#ifdef MICM_USE_CUDA
  template<micm::Index L>
  using CudaRosen = micm::CudaSolverBuilderInPlace<micm::RosenbrockSolverParameters, L>;
#endif
#ifdef MICM_USE_KOKKOS
  template<micm::Index L>
  using KokkosRosen = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters, L>;
#endif

  /// @brief Register the four CPU LU combinations for one matrix ordering.
  template<class Mechanism, class DM, class SM>
  void RegisterCpu(Registry& registry, const std::string& matrix)
  {
    const std::string mechanism{ Mechanism::kName };
    registry[Key(mechanism, "cpu", matrix, "in-place", "mozart")] =
        &RunCase<Mechanism, CpuRosenInPlace<DM, SM, MozLUInPlace<SM>>>;
    registry[Key(mechanism, "cpu", matrix, "in-place", "doolittle")] =
        &RunCase<Mechanism, CpuRosenInPlace<DM, SM, DooLUInPlace<SM>>>;
    registry[Key(mechanism, "cpu", matrix, "separate", "mozart")] = &RunCase<Mechanism, CpuRosen<DM, SM, MozLU<SM>>>;
    registry[Key(mechanism, "cpu", matrix, "separate", "doolittle")] = &RunCase<Mechanism, CpuRosen<DM, SM, DooLU<SM>>>;
  }

  /// @brief Register every configuration that a single vector width supports.
  ///        The CUDA backend only supports the in-place Mozart LU.
  ///        Kokkos backend supports all LU types, but we just include Mozart LU to cut down on build times
  template<class Mechanism, micm::Index L>
  void RegisterVector(Registry& registry)
  {
    const std::string matrix = "vector" + std::to_string(L);
    RegisterCpu<Mechanism, VectorDense<L>, VectorSparse<L>>(registry, matrix);
#ifdef MICM_USE_CUDA
    registry[Key(std::string{ Mechanism::kName }, "gpu", matrix, "in-place", "mozart")] = &RunCase<Mechanism, CudaRosen<L>>;
#endif
#ifdef MICM_USE_KOKKOS
    registry[Key(std::string{ Mechanism::kName }, "kokkos", matrix, "in-place", "mozart")] =
        &RunCase<Mechanism, KokkosRosen<L>>;
#endif
  }

  /// @brief Register one mechanism against every solver configuration.
  template<class Mechanism, micm::Index... Ls>
  void Register(Registry& registry, std::integer_sequence<micm::Index, Ls...>)
  {
    RegisterCpu<Mechanism, StandardDense, StandardSparse>(registry, "standard");
    (RegisterVector<Mechanism, Ls>(registry), ...);
  }

  /// @brief The vector widths the benchmark reports on.
  using VectorWidths = std::integer_sequence<micm::Index, 1, 2, 4, 8, 128>;

}  // namespace bench
