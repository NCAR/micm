#include <chrono>
#include <iomanip>
#include <iostream>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

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
      auto result = solver.template Solve<true>(time_step - elapsed_solve_time, state);
      auto end = std::chrono::high_resolution_clock::now();

      total_solve_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
      total_stats.function_calls += result.stats_.function_calls;
      total_stats.jacobian_updates += result.stats_.jacobian_updates;
      total_stats.number_of_steps += result.stats_.number_of_steps;
      total_stats.accepted += result.stats_.accepted;
      total_stats.rejected += result.stats_.rejected;
      total_stats.decompositions += result.stats_.decompositions;
      total_stats.solves += result.stats_.solves;
      total_stats.singular += result.stats_.singular;
      total_stats.total_forcing_time += result.stats_.total_forcing_time;
      total_stats.total_jacobian_time += result.stats_.total_jacobian_time;
      total_stats.total_linear_factor_time += result.stats_.total_linear_factor_time;
      total_stats.total_linear_solve_time += result.stats_.total_linear_solve_time;

      elapsed_solve_time = result.final_time_;
      state.variables_ = result.result_;
    }

    state.PrintState(time_step * (i + 1));
  }
  std::cout << "Total solve time: " << total_solve_time.count() << " nanoseconds" << std::endl;
  std::cout << "accepted: " << total_stats.accepted << std::endl;
  std::cout << "function_calls: " << total_stats.function_calls << std::endl;
  std::cout << "jacobian_updates: " << total_stats.jacobian_updates << std::endl;
  std::cout << "number_of_steps: " << total_stats.number_of_steps << std::endl;
  std::cout << "accepted: " << total_stats.accepted << std::endl;
  std::cout << "rejected: " << total_stats.rejected << std::endl;
  std::cout << "decompositions: " << total_stats.decompositions << std::endl;
  std::cout << "solves: " << total_stats.solves << std::endl;
  std::cout << "singular: " << total_stats.singular << std::endl;
  std::cout << "total_forcing_time: " << total_stats.total_forcing_time.count() << " nanoseconds" << std::endl;
  std::cout << "total_jacobian_time: " << total_stats.total_jacobian_time.count() << " nanoseconds" << std::endl;
  std::cout << "total_linear_factor_time: " << total_stats.total_linear_factor_time.count() << " nanoseconds" << std::endl;
  std::cout << "total_linear_solve_time: " << total_stats.total_linear_solve_time.count() << " nanoseconds" << std::endl;
}

int main()
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species("C");

  Phase gas_phase{ std::vector<Species>{ a, b, c } };

  Process r1 = Process::create()
                   .reactants({ a })
                   .products({ yields(b, 1) })
                   .rate_constant(UserDefinedRateConstant({ .label_ = "r1" }))
                   .phase(gas_phase);

  Process r2 = Process::create()
                   .reactants({ b, b })
                   .products({ yields(b, 1), yields(c, 1) })
                   .rate_constant(UserDefinedRateConstant({ .label_ = "r2" }))
                   .phase(gas_phase);

  Process r3 = Process::create()
                   .reactants({ b, c })
                   .products({ yields(a, 1), yields(c, 1) })
                   .rate_constant(UserDefinedRateConstant({ .label_ = "r3" }))
                   .phase(gas_phase);

  auto system = System(SystemParameters{ .gas_phase_ = gas_phase });
  auto reactions = std::vector<Process>{ r1, r2, r3 };

  RosenbrockSolver<> two_stage{ system, reactions, RosenbrockSolverParameters::two_stage_rosenbrock_parameters() };

  RosenbrockSolver<> three_stage{ system, reactions, RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };

  RosenbrockSolver<> four_stage{ system, reactions, RosenbrockSolverParameters::four_stage_rosenbrock_parameters() };

  RosenbrockSolver<> four_stage_da{ system,
                                    reactions,
                                    RosenbrockSolverParameters::four_stage_differential_algebraic_rosenbrock_parameters() };

  RosenbrockSolver<> six_stage_da{ system,
                                   reactions,
                                   RosenbrockSolverParameters::six_stage_differential_algebraic_rosenbrock_parameters() };

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