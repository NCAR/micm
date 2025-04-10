#include <micm/JIT.hpp>
#include <chrono>
#include <iostream>

// Use our namespace so that this example is easier to read
using namespace micm;

constexpr size_t n_grid_cells = 1;

// partial template specializations
template<class T>
using GroupVectorMatrix = micm::VectorMatrix<T, n_grid_cells>;
using GroupSparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<n_grid_cells>>;

auto run_solver(auto& solver)
{
  SolverStats total_stats;
  State state = solver.GetState();

  state.variables_ = 1;

  for (int i = 0; i < n_grid_cells; ++i)
  {
    state.conditions_[i].temperature_ = 287.45;  // K
    state.conditions_[i].pressure_ = 101319.9;   // Pa
    state.conditions_[i].air_density_ = 1e6;     // mol m-3
  }
  auto foo = Species("Foo");
  std::vector<double> foo_conc(n_grid_cells, 1.0);
  state.SetConcentration(foo, foo_conc);

  // choose a timestep and print the initial state
  double time_step = 500;  // s

  auto total_solve_time = std::chrono::nanoseconds::zero();

  // solve for ten iterations
  for (int i = 0; i < 10; ++i)
  {
    double elapsed_solve_time = 0;

    while (elapsed_solve_time < time_step)
    {
      auto start = std::chrono::high_resolution_clock::now();
      solver.CalculateRateConstants(state);
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      auto end = std::chrono::high_resolution_clock::now();
      total_solve_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
      elapsed_solve_time = result.final_time_;

      total_stats.function_calls_ += result.stats_.function_calls_;
      total_stats.jacobian_updates_ += result.stats_.jacobian_updates_;
      total_stats.number_of_steps_ += result.stats_.number_of_steps_;
      total_stats.accepted_ += result.stats_.accepted_;
      total_stats.rejected_ += result.stats_.rejected_;
      total_stats.decompositions_ += result.stats_.decompositions_;
      total_stats.solves_ += result.stats_.solves_;
    }
  }

  return std::make_tuple(state, total_stats, total_solve_time);
}

int main(const int argc, const char* argv[])
{
  auto foo = Species{ "Foo" };
  auto bar = Species{ "Bar" };
  auto baz = Species{ "Baz" };

  Phase gas_phase{ std::vector<Species>{ foo, bar, baz } };

  System chemical_system{ SystemParameters{ .gas_phase_ = gas_phase } };

  Process r1 = Process::Create()
                   .SetReactants({ foo })
                   .SetProducts({ Yield(bar, 0.8), Yield(baz, 0.2) })
                   .SetRateConstant(ArrheniusRateConstant({ .A_ = 1.0e-3 }))
                   .SetPhase(gas_phase);

  Process r2 = Process::Create()
                   .SetReactants({ foo, bar })
                   .SetProducts({ Yield(baz, 1) })
                   .SetRateConstant(ArrheniusRateConstant({ .A_ = 1.0e-5, .C_ = 110.0 }))
                   .SetPhase(gas_phase);

  std::vector<Process> reactions{ r1, r2 };

  auto solver_parameters = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();

  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(solver_parameters)
                    .SetSystem(chemical_system)
                    .SetReactions(reactions)
                    .SetNumberOfGridCells(n_grid_cells)
                    .Build();
  
  auto start = std::chrono::high_resolution_clock::now();
  auto jit_solver = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, n_grid_cells>(solver_parameters)
                        .SetSystem(chemical_system)
                        .SetReactions(reactions)
                        .SetNumberOfGridCells(n_grid_cells)
                        .Build();
  auto end = std::chrono::high_resolution_clock::now();
  auto jit_compile_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

  std::cout << "Jit compile time: " << jit_compile_time.count() << " nanoseconds" << std::endl;

  auto result_tuple = run_solver(solver);
  auto jit_result_tuple = run_solver(jit_solver);

  // Rerun for more fair comparison after assumed improvements to
  // branch-prediction during state update
  result_tuple = run_solver(solver);
  jit_result_tuple = run_solver(jit_solver);

  std::cout << "Standard solve time: " << std::get<2>(result_tuple).count() << " nanoseconds" << std::endl;
  std::cout << "JIT solve time: " << std::get<2>(jit_result_tuple).count() << " nanoseconds" << std::endl;

  auto result_stats = std::get<1>(result_tuple);
  std::cout << "Standard solve stats: " << std::endl;
  std::cout << "\taccepted: " << result_stats.accepted_ << std::endl;
  std::cout << "\tfunction_calls: " << result_stats.function_calls_ << std::endl;
  std::cout << "\tjacobian_updates: " << result_stats.jacobian_updates_ << std::endl;
  std::cout << "\tnumber_of_steps: " << result_stats.number_of_steps_ << std::endl;
  std::cout << "\taccepted: " << result_stats.accepted_ << std::endl;
  std::cout << "\trejected: " << result_stats.rejected_ << std::endl;
  std::cout << "\tdecompositions: " << result_stats.decompositions_ << std::endl;
  std::cout << "\tsolves: " << result_stats.solves_ << std::endl;

  auto jit_result_stats = std::get<1>(jit_result_tuple);
  std::cout << "JIT solve stats: " << std::endl;
  std::cout << "\taccepted: " << jit_result_stats.accepted_ << std::endl;
  std::cout << "\tfunction_calls: " << jit_result_stats.function_calls_ << std::endl;
  std::cout << "\tjacobian_updates: " << jit_result_stats.jacobian_updates_ << std::endl;
  std::cout << "\tnumber_of_steps: " << jit_result_stats.number_of_steps_ << std::endl;
  std::cout << "\taccepted: " << jit_result_stats.accepted_ << std::endl;
  std::cout << "\trejected: " << jit_result_stats.rejected_ << std::endl;
  std::cout << "\tdecompositions: " << jit_result_stats.decompositions_ << std::endl;
  std::cout << "\tsolves: " << jit_result_stats.solves_ << std::endl;

  auto result = std::get<0>(result_tuple);
  auto jit_result = std::get<0>(jit_result_tuple);

  for (auto& species : result.variable_names_)
  {
    for (int i = 0; i < n_grid_cells; ++i)
    {
      double a = result.variables_[i][result.variable_map_[species]];
      double b = jit_result.variables_[i][jit_result.variable_map_[species]];
      if (std::abs(a - b) > 1.0e-5 * (std::abs(a) + std::abs(b)) / 2.0 + 1.0e-30)
      {
        std::cout << species << " does not match final concentration" << std::endl;
      }
    }
  }
  return 0;
}