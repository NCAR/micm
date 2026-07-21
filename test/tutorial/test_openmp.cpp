#include <micm/CPU.hpp>
#include <micm/util/types.hpp>

#include <omp.h>

#include <iomanip>

using namespace micm;

void print_header()
{
  std::cout << std::setw(10) << "A"
            << "," << std::setw(10) << "B"
            << "," << std::setw(10) << "C" << std::endl;
}

void print_results(std::vector<Real> results)
{
  std::ios oldState(nullptr);
  oldState.copyfmt(std::cout);

  std::cout << std::scientific << std::setprecision(2) << std::setw(10) << results[0] << "," << std::setw(10) << results[1]
            << "," << std::setw(10) << results[2] << std::endl;

  std::cout.copyfmt(oldState);
}

std::vector<Real> run_solver_on_thread_with_own_state(auto& solver, auto& state)
{
  std::cout << "Running solver on thread " << omp_get_thread_num() << std::endl;

  // mol m-3
  state.variables_[0] = { 1, 0, 0 };

  Real k1 = 0.04;
  Real k2 = 3e7;
  Real k3 = 1e4;
  state.SetCustomRateParameter("r1", k1);
  state.SetCustomRateParameter("r2", k2);
  state.SetCustomRateParameter("r3", k3);

  Real temperature = 272.5;  // [K]
  Real pressure = 101253.3;  // [Pa]
  Real air_density = 1e6;    // [mol m-3]

  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  Real time_step = 200;  // s

  for (Index i = 0; i < 10; ++i)
  {
    Real elapsed_solve_time = 0;
    solver.UpdateStateParameters(state);

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.stats_.final_time_;;
    }
  }

  return state.variables_.AsVector();
}

int main()
{
  constexpr Index n_threads = 3;

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ a, b, c } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ micm::StoichSpecies(b, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r1" })
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b, b })
                         .SetProducts({ micm::StoichSpecies(b, 1), micm::StoichSpecies(c, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r2" })
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b, c })
                         .SetProducts({ micm::StoichSpecies(a, 1), micm::StoichSpecies(c, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "r3" })
                         .SetPhase(gas_phase)
                         .Build();

  auto reactions = std::vector<micm::Process>{ r1, r2, r3 };
  auto chemical_system = micm::System(gas_phase );

  std::vector<std::vector<Real>> results(n_threads);

  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
                    micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(chemical_system)
                    .SetReactions(reactions)
                    .Build();

#pragma omp parallel num_threads(n_threads)
  {
    // each thread should use its own state
    auto state = solver.GetState();
    std::vector<Real> result = run_solver_on_thread_with_own_state(solver, state);
    results[omp_get_thread_num()] = result;
#pragma omp barrier
  }

  std::cout << "Thread 1" << std::endl;
  print_header();
  print_results(results[0]);

  std::cout << "Thread 2" << std::endl;
  print_header();
  print_results(results[1]);

  std::cout << "Thread 3" << std::endl;
  print_header();
  print_results(results[2]);
}