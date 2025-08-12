#include <micm/CPU.hpp>

#include <chrono>
#include <iomanip>
#include <iostream>

// Use our namespace so that this example is easier to read
using namespace micm;

void test_solver_type(auto& solver)
{
  auto state = solver.GetState();

  // mol m-3
  state.variables_[0] = { 1, 0, 0 };

  double k1 = 0.04;
  double k2 = 3e7;
  double k3 = 1e4;
  state.SetCustomRateParameter("r1", k1);
  state.SetCustomRateParameter("r2", k2);
  state.SetCustomRateParameter("r3", k3);

  double temperature = 272.5;  // [K]
  double pressure = 101253.3;  // [Pa]
  double air_density = 1e6;    // [mol m-3]

  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.conditions_[0].air_density_ = air_density;

  // choose a timestep and print the initial state
  double time_step = 200;  // s

  state.PrintHeader();
  state.PrintState(0);

  SolverStats total_stats;
  std::chrono::duration<double, std::nano> total_solve_time = std::chrono::nanoseconds::zero();

  // solve for ten iterations
  for (int i = 0; i < 10; ++i)
  {
    // Depending on how stiff the system is
    // the solver integration step may not be able to solve for the full time step
    // so we need to track how much time the solver was able to integrate for and continue
    // solving until we finish
    double elapsed_solve_time = 0;

    while (elapsed_solve_time < time_step)
    {
      auto start = std::chrono::high_resolution_clock::now();
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      auto end = std::chrono::high_resolution_clock::now();

      total_solve_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
      total_stats.function_calls_ += result.stats_.function_calls_;
      total_stats.jacobian_updates_ += result.stats_.jacobian_updates_;
      total_stats.number_of_steps_ += result.stats_.number_of_steps_;
      total_stats.accepted_ += result.stats_.accepted_;
      total_stats.rejected_ += result.stats_.rejected_;
      total_stats.decompositions_ += result.stats_.decompositions_;
      total_stats.solves_ += result.stats_.solves_;

      elapsed_solve_time += result.final_time_;
    }

    state.PrintState(time_step * (i + 1));
  }
  std::cout << "Total solve time: " << total_solve_time.count() << " nanoseconds" << std::endl;
  std::cout << "accepted: " << total_stats.accepted_ << std::endl;
  std::cout << "function_calls: " << total_stats.function_calls_ << std::endl;
  std::cout << "jacobian_updates: " << total_stats.jacobian_updates_ << std::endl;
  std::cout << "number_of_steps: " << total_stats.number_of_steps_ << std::endl;
  std::cout << "accepted: " << total_stats.accepted_ << std::endl;
  std::cout << "rejected: " << total_stats.rejected_ << std::endl;
  std::cout << "decompositions: " << total_stats.decompositions_ << std::endl;
  std::cout << "solves: " << total_stats.solves_ << std::endl;
}

int main()
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species("C");

  Phase gas_phase{ std::vector<Species>{ a, b, c } };

  Process r1 = ChemicalReactionBuilder()
                   .SetReactants({ a })
                   .SetProducts({ Yield(b, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "r1" }))
                   .SetPhaseName("gas")
                   .Build();

  Process r2 = ChemicalReactionBuilder()
                   .SetReactants({ b, b })
                   .SetProducts({ Yield(b, 1), Yield(c, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "r2" }))
                   .SetPhaseName("gas")
                   .Build();

  Process r3 = ChemicalReactionBuilder()
                   .SetReactants({ b, c })
                   .SetProducts({ Yield(a, 1), Yield(c, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "r3" }))
                   .SetPhaseName("gas")
                   .Build();

  auto system = System(SystemParameters{ .gas_phase_ = gas_phase });
  auto reactions = std::vector<Process>{ r1, r2, r3 };

  auto two_stage = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
                       micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters())
                       .SetSystem(system)
                       .SetReactions(reactions)
                       .Build();

  auto three_stage = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
                         micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                         .SetSystem(system)
                         .SetReactions(reactions)
                         .Build();

  auto four_stage = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
                        micm::RosenbrockSolverParameters::FourStageRosenbrockParameters())
                        .SetSystem(system)
                        .SetReactions(reactions)
                        .Build();

  auto four_stage_da = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
                           micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
                           .SetSystem(system)
                           .SetReactions(reactions)
                           .Build();

  auto six_stage_da = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
                          micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters())
                          .SetSystem(system)
                          .SetReactions(reactions)
                          .Build();

  std::cout << "Two stages: " << std::endl;
  test_solver_type(two_stage);

  std::cout << std::endl << "Three stages: " << std::endl;
  test_solver_type(three_stage);

  std::cout << std::endl << "Four stages: " << std::endl;
  test_solver_type(four_stage);

  std::cout << std::endl << "Four stages differential algebraic: " << std::endl;
  test_solver_type(four_stage_da);

  std::cout << std::endl << "Six stages differential algebraic: " << std::endl;
  test_solver_type(six_stage_da);
}