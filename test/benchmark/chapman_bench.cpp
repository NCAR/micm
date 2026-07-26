// Standalone Chapman benchmark for measuring solver performance.
// Runs the same 7-reaction Chapman mechanism as
// test/integration/test_chapman_integration.cpp, but scaled to many grid cells
// so per-Solve() work dominates over per-call overhead.
//
// Usage: chapman_bench [num_cells] [num_steps] [dt_seconds] [matrix_kind]
//   matrix_kind: "standard" (default) or "vector1|2|4|8|128"
//
// Prints a single line with configuration and elapsed wall time (ms).

#include <micm/CPU.hpp>

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

#define SOLVER_BUILDER micm::CpuSolverBuilder

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
  void InitState(typename Solver::StatePolicyType& state, std::size_t num_cells)
  {
    // Initial concentrations: matches the integration test's cell 0 for consistency.
    std::vector<double> concentrations{ 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3 };
    std::vector<double> photo_rates{ 0.1, 0.2, 0.3 };
    for (std::size_t c = 0; c < num_cells; ++c)
    {
      state.variables_[c] = concentrations;
      state.custom_rate_parameters_[c] = photo_rates;
      state.conditions_[c].temperature_ = 273.15 + 25.0;
      state.conditions_[c].pressure_ = 101325.0;
    }
  }

  template<class Solver>
  double RunBench(Solver& solver, std::size_t num_cells, std::size_t num_steps, double dt)
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
    for (std::size_t s = 0; s < num_steps; ++s)
    {
      [[maybe_unused]] auto r = solver.Solve(dt, state);
    }
    auto t1 = std::chrono::steady_clock::now();
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
  }
}  // namespace

int main(int argc, char** argv)
{
  std::size_t num_cells = (argc > 1) ? std::stoul(argv[1]) : 10000;
  std::size_t num_steps = (argc > 2) ? std::stoul(argv[2]) : 100;
  double dt = (argc > 3) ? std::stod(argv[3]) : 30.0;
  std::string kind = (argc > 4) ? argv[4] : "standard";

  double elapsed_ms = 0.0;
  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();

  if (kind == "standard")
  {
    auto solver =
        BuildChapmanSolver(SOLVER_BUILDER<micm::RosenbrockSolverParameters>(options));
    elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
  }
  else if (kind == "vector1")
  {
    auto solver = BuildChapmanSolver(SOLVER_BUILDER<
                                     micm::RosenbrockSolverParameters,
                                     micm::VectorMatrix<double, 1>,
                                     micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>>(options));
    elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
  }
  else if (kind == "vector2")
  {
    auto solver = BuildChapmanSolver(SOLVER_BUILDER<
                                     micm::RosenbrockSolverParameters,
                                     micm::VectorMatrix<double, 2>,
                                     micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>>(options));
    elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
  }
  else if (kind == "vector4")
  {
    auto solver = BuildChapmanSolver(SOLVER_BUILDER<
                                     micm::RosenbrockSolverParameters,
                                     micm::VectorMatrix<double, 4>,
                                     micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<4>>>(options));
    elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
  }
  else if (kind == "vector8")
  {
    auto solver = BuildChapmanSolver(SOLVER_BUILDER<
                                     micm::RosenbrockSolverParameters,
                                     micm::VectorMatrix<double, 8>,
                                     micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<8>>>(options));
    elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
  }
  else if (kind == "vector128")
  {
    auto solver = BuildChapmanSolver(SOLVER_BUILDER<
                                     micm::RosenbrockSolverParameters,
                                     micm::VectorMatrix<double, 128>,
                                     micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<128>>>(options));
    elapsed_ms = RunBench(solver, num_cells, num_steps, dt);
  }
  else
  {
    std::cerr << "Unknown matrix kind: " << kind
              << " (expected standard|vector1|vector2|vector4|vector8|vector128)\n";
    return 1;
  }

  std::cout << "kind=" << kind << " cells=" << num_cells << " steps=" << num_steps
            << " dt=" << dt << " elapsed_ms=" << elapsed_ms
            << " ms_per_step=" << elapsed_ms / static_cast<double>(num_steps) << "\n";
  return 0;
}
