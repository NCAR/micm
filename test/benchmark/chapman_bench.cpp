// Standalone Chapman benchmark for measuring solver performance.
// Runs the same 7-reaction Chapman mechanism as
// test/integration/test_chapman_integration.cpp, but scaled to many grid cells
// so per-Solve() work dominates over per-call overhead.
//
// Usage: chapman_bench [num_cells] [num_steps] [dt_seconds] [backend] [matrix] [lu_type] [lu_algorithm]
//   backend: "cpu" (default) or "gpu"
//   matrix: "standard" (default) or "vector1|2|4|8|128"
//   lu_type: "in-place" (default) or "separate"
//   lu_algorithm: "mozart" (default) or "doolittle"
//
//   Note that "gpu" backend must use "mozart"/"in-place" LU and "vector1|2|4|8|128" matrix
//
// Prints a single line with configuration and elapsed wall time (ms).

#include <micm/CPU.hpp>
#ifdef MICM_USE_CUDA
#include <micm/GPU.hpp>
#include <micm/cuda/util/cuda_util.cuh>
#endif
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// Callgrind client-request macros. When compiled without valgrind-devel
// available (or when the binary is run outside valgrind), these are no-ops,
// so the benchmark works the same under `./chapman_bench` and under
// `valgrind --tool=callgrind ./chapman_bench`.
//
// Under callgrind the profile script invokes us with
//     --collect-atstart=no --instr-atstart=no
// so instrumentation is completely disabled during solver construction, state
// initialization, rate-constant setup, and the warm-up call. We toggle it on
// only around the timed Solve() loop, so the reported instruction count
// covers exactly Rosenbrock Solve() and everything it calls.
#if __has_include(<valgrind/callgrind.h>)
  #include <valgrind/callgrind.h>
#else
  #define CALLGRIND_START_INSTRUMENTATION do {} while (0)
  #define CALLGRIND_STOP_INSTRUMENTATION  do {} while (0)
  #define CALLGRIND_TOGGLE_COLLECT        do {} while (0)
  #define CALLGRIND_ZERO_STATS            do {} while (0)
#endif

namespace
{
  template<class Builder>
  auto BuildChapmanSolver(Builder builder)
  {
    auto o = micm::Species("O");
    auto o1d = micm::Species("O1D");
    auto o2 = micm::Species("O2");
    auto o3 = micm::Species("O3");
    auto m = micm::Species("M");
    auto ar = micm::Species("Ar");
    auto n2 = micm::Species("N2");
    auto h2o = micm::Species("H2O");
    auto co2 = micm::Species("CO2");

    micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ o, o1d, o2, o3, m, ar, n2, h2o, co2 } };

    micm::Process r1 = micm::ChemicalReactionBuilder()
                           .SetReactants({ o1d, n2 })
                           .SetProducts({ micm::StoichSpecies(o, 1), micm::StoichSpecies(n2, 1) })
                           .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 })
                           .SetPhase(gas_phase)
                           .Build();

    micm::Process r2 = micm::ChemicalReactionBuilder()
                           .SetReactants({ o1d, o2 })
                           .SetProducts({ micm::StoichSpecies(o, 1), micm::StoichSpecies(o2, 1) })
                           .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 })
                           .SetPhase(gas_phase)
                           .Build();

    micm::Process r3 = micm::ChemicalReactionBuilder()
                           .SetReactants({ o, o3 })
                           .SetProducts({ micm::StoichSpecies(o2, 2) })
                           .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 8e-12, .B_ = 0, .C_ = -2060 })
                           .SetPhase(gas_phase)
                           .Build();

    micm::Process r4 = micm::ChemicalReactionBuilder()
                           .SetReactants({ o, o2, m })
                           .SetProducts({ micm::StoichSpecies(o3, 1), micm::StoichSpecies(m, 1) })
                           .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6.0e-34, .B_ = 0, .C_ = 2.4 })
                           .SetPhase(gas_phase)
                           .Build();

    micm::Process photo_1 = micm::ChemicalReactionBuilder()
                                .SetReactants({ o2 })
                                .SetProducts({ micm::StoichSpecies(o, 2) })
                                .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jO2" })
                                .SetPhase(gas_phase)
                                .Build();

    micm::Process photo_2 = micm::ChemicalReactionBuilder()
                                .SetReactants({ o3 })
                                .SetProducts({ micm::StoichSpecies(o1d, 1), micm::StoichSpecies(o2, 1) })
                                .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jO3a" })
                                .SetPhase(gas_phase)
                                .Build();

    micm::Process photo_3 = micm::ChemicalReactionBuilder()
                                .SetReactants({ o3 })
                                .SetProducts({ micm::StoichSpecies(o, 1), micm::StoichSpecies(o2, 1) })
                                .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "jO3b" })
                                .SetPhase(gas_phase)
                                .Build();

    return builder.SetSystem(micm::System(gas_phase))
        .SetReactions({ r1, r2, r3, r4, photo_1, photo_2, photo_3 })
        .SetIgnoreUnusedSpecies(true)
        .Build();
  }

  template<class Solver>
  void InitState(typename Solver::StatePolicyType& state, micm::Index num_cells)
  {
    // Initial concentrations: matches the integration test's cell 0 for consistency.
    std::vector<micm::Real> concentrations{ 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3 };
    std::vector<micm::Real> photo_rates{ 0.1, 0.2, 0.3 };
    for (micm::Index c = 0; c < num_cells; ++c)
    {
      state.variables_[c] = concentrations;
      state.custom_rate_parameters_[c] = photo_rates;
      state.conditions_[c].temperature_ = 273.15 + 25.0;
      state.conditions_[c].pressure_ = 101325.0;
    }
  }

  template<class Solver>
  double RunBench(Solver& solver, micm::Index num_cells, micm::Index num_steps, micm::Real dt)
  {
    auto state = solver.GetState(num_cells);
    InitState<Solver>(state, num_cells);
    solver.UpdateStateParameters(state);

    // Warm-up (first call primes caches / any lazy init).
    [[maybe_unused]] auto warmup = solver.Solve(dt, state);

    // Turn callgrind instrumentation on just for the timed loop. When not
    // running under callgrind these macros are no-ops.
    CALLGRIND_ZERO_STATS;
    CALLGRIND_START_INSTRUMENTATION;
    CALLGRIND_TOGGLE_COLLECT;
    auto t0 = std::chrono::steady_clock::now();
    for (micm::Index s = 0; s < num_steps; ++s)
    {
      [[maybe_unused]] auto r = solver.Solve(dt, state);
    }
    auto t1 = std::chrono::steady_clock::now();
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
  }

  using StandardDense = micm::Matrix<micm::Real>;
  template<micm::Index L>
  using VectorDense = micm::VectorMatrix<micm::Real, L>;
  using StandardSparse = micm::SparseMatrix<micm::Real, micm::SparseMatrixStandardOrdering>;
  template<micm::Index L>
  using VectorSparse = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<L>>;

  template<class SparseMatrixPolicy>
  using DooLU = micm::LuDecompositionDoolittle<SparseMatrixPolicy>;
  template<class SparseMatrixPolicy>
  using DooLUInPlace = micm::LuDecompositionDoolittleInPlace<SparseMatrixPolicy>;
  template<class SparseMatrixPolicy>
  using MozLU = micm::LuDecompositionMozart<SparseMatrixPolicy>;
  template<class SparseMatrixPolicy>
  using MozLUInPlace = micm::LuDecompositionMozartInPlace<SparseMatrixPolicy>;

  template<class DM, class SM, class LU>
  using CpuRosen = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters, DM, SM, LU>;
  template<class DM, class SM, class LU>
  using CpuRosenInPlace = micm::CpuSolverBuilderInPlace<micm::RosenbrockSolverParameters, DM, SM, LU>;
#ifdef MICM_USE_CUDA
  template<micm::Index L>
  using CudaRosen = micm::CudaSolverBuilderInPlace<micm::RosenbrockSolverParameters, L>;
#endif
}  // namespace

int main(int argc, char** argv)
{
  micm::Index num_cells = (argc > 1) ? std::stoul(argv[1]) : 10000;
  micm::Index num_steps = (argc > 2) ? std::stoul(argv[2]) : 100;
  micm::Real dt = (argc > 3) ? static_cast<micm::Real>(std::stod(argv[3])) : static_cast<micm::Real>(30.0);
  std::string backend = (argc > 4) ? argv[4] : "cpu"; // "cpu", "gpu"
  std::string matrix_type = (argc > 5) ? argv[5] : "standard"; // "standard", "vector1", "vector2", "vector4", "vector8", "vector128"
  std::string lu_matrix_type = (argc > 6) ? argv[6] : "in-place"; // "in-place", "separate"
  std::string lu_type = (argc > 7) ? argv[7] : "mozart"; // "mozart", "doolittle"

  double elapsed_ms = -1.0;
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();

  if (backend == "cpu")
  {
    if (matrix_type == "standard")
    {
      if (lu_matrix_type == "in-place")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<StandardDense, StandardSparse, MozLUInPlace<StandardSparse>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<StandardDense, StandardSparse, DooLUInPlace<StandardSparse>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
      else if (lu_matrix_type == "separate")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosen<StandardDense, StandardSparse, MozLU<StandardSparse>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosen<StandardDense, StandardSparse, DooLU<StandardSparse>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
    }
    else if (matrix_type == "vector1")
    {
      if (lu_matrix_type == "in-place")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<VectorDense<1>, VectorSparse<1>, MozLUInPlace<VectorSparse<1>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<VectorDense<1>, VectorSparse<1>, DooLUInPlace<VectorSparse<1>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
      else if (lu_matrix_type == "separate")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosen<VectorDense<1>, VectorSparse<1>, MozLU<VectorSparse<1>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosen<VectorDense<1>, VectorSparse<1>, DooLU<VectorSparse<1>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
    }
    else if (matrix_type == "vector2")
    {
      if (lu_matrix_type == "in-place")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<VectorDense<2>, VectorSparse<2>, MozLUInPlace<VectorSparse<2>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<VectorDense<2>, VectorSparse<2>, DooLUInPlace<VectorSparse<2>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
      else if (lu_matrix_type == "separate")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosen<VectorDense<2>, VectorSparse<2>, MozLU<VectorSparse<2>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosen<VectorDense<2>, VectorSparse<2>, DooLU<VectorSparse<2>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
    }
    else if (matrix_type == "vector4")
    {
      if (lu_matrix_type == "in-place")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<VectorDense<4>, VectorSparse<4>, MozLUInPlace<VectorSparse<4>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<VectorDense<4>, VectorSparse<4>, DooLUInPlace<VectorSparse<4>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
      else if (lu_matrix_type == "separate")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosen<VectorDense<4>, VectorSparse<4>, MozLU<VectorSparse<4>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosen<VectorDense<4>, VectorSparse<4>, DooLU<VectorSparse<4>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
    }
    else if (matrix_type == "vector8")
    {
      if (lu_matrix_type == "in-place")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<VectorDense<8>, VectorSparse<8>, MozLUInPlace<VectorSparse<8>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<VectorDense<8>, VectorSparse<8>, DooLUInPlace<VectorSparse<8>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
      else if (lu_matrix_type == "separate")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosen<VectorDense<8>, VectorSparse<8>, MozLU<VectorSparse<8>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosen<VectorDense<8>, VectorSparse<8>, DooLU<VectorSparse<8>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
    }
    else if (matrix_type == "vector128")
    {
      if (lu_matrix_type == "in-place")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<VectorDense<128>, VectorSparse<128>, MozLUInPlace<VectorSparse<128>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosenInPlace<VectorDense<128>, VectorSparse<128>, DooLUInPlace<VectorSparse<128>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
      else if (lu_matrix_type == "separate")
      {
        if (lu_type == "mozart")
        {
          auto solver = BuildChapmanSolver(CpuRosen<VectorDense<128>, VectorSparse<128>, MozLU<VectorSparse<128>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
        else if (lu_type == "doolittle")
        {
          auto solver = BuildChapmanSolver(CpuRosen<VectorDense<128>, VectorSparse<128>, DooLU<VectorSparse<128>>>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
        }
      }
    }
  }
#ifdef MICM_USE_CUDA
  else if (backend == "gpu" && lu_matrix_type == "in-place" && lu_type == "mozart")
  {
    if (matrix_type == "vector1")
    {
      auto solver = BuildChapmanSolver(CudaRosen<1>(options));
      elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
    }
    else if (matrix_type == "vector2")
    {
      auto solver = BuildChapmanSolver(CudaRosen<2>(options));
      elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
    }
    else if (matrix_type == "vector4")
    {
          auto solver = BuildChapmanSolver(CudaRosen<4>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
    }
    else if (matrix_type == "vector8")
    {
          auto solver = BuildChapmanSolver(CudaRosen<8>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
    }
    else if (matrix_type == "vector128")
    {
          auto solver = BuildChapmanSolver(CudaRosen<128>(options));
          elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
    }
  }
  micm::cuda::CudaStreamSingleton::GetInstance().CleanUp();
#endif
  if (elapsed_ms < 0.0)
  {
    std::cout << "Invalid option combination: backend='" << backend << "'; matrix_type='" << matrix_type
            << "'; lu_matrix_type='" << lu_matrix_type << "'; lu_type='" << lu_type << "'" << std::endl;
    return 1;
  }

  std::cout << "backend=" << backend << " matrix_type=" << matrix_type << " lu=" << lu_type << "/"
            << lu_matrix_type << " cells=" << num_cells << " steps=" << num_steps
            << " dt=" << dt << " elapsed_ms=" << elapsed_ms
            << " ms_per_step=" << elapsed_ms / static_cast<double>(num_steps) << std::endl;
  return 0;
}
