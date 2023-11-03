#include <chrono>
#include <iostream>
#include <micm/configure/solver_config.hpp>
#include <micm/jit/jit_compiler.hpp>
#include <micm/solver/jit_rosenbrock.hpp>
#include <micm/solver/rosenbrock.hpp>

// Use our namespace so that this example is easier to read
using namespace micm;

constexpr size_t n_grid_cells = 1;

// partial template specializations
template<class T>
using GroupVectorMatrix = micm::VectorMatrix<T, n_grid_cells>;
template<class T>
using GroupSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<n_grid_cells>>;

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

  state.SetCustomRateParameter("EMIS.NO", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.44e-10));
  state.SetCustomRateParameter("EMIS.NO2", std::vector<double>(solver.parameters_.number_of_grid_cells_, 7.56e-12));
  state.SetCustomRateParameter("EMIS.CO", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.96e-09));
  state.SetCustomRateParameter("EMIS.SO2", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.06e-09));
  state.SetCustomRateParameter("EMIS.FORM", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.02e-11));
  state.SetCustomRateParameter("EMIS.MEOH", std::vector<double>(solver.parameters_.number_of_grid_cells_, 5.92e-13));
  state.SetCustomRateParameter("EMIS.ALD2", std::vector<double>(solver.parameters_.number_of_grid_cells_, 4.25e-12));
  state.SetCustomRateParameter("EMIS.PAR", std::vector<double>(solver.parameters_.number_of_grid_cells_, 4.27e-10));
  state.SetCustomRateParameter("EMIS.ETH", std::vector<double>(solver.parameters_.number_of_grid_cells_, 4.62e-11));
  state.SetCustomRateParameter("EMIS.OLE", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.49e-11));
  state.SetCustomRateParameter("EMIS.IOLE", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.49e-11));
  state.SetCustomRateParameter("EMIS.TOL", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.53e-11));
  state.SetCustomRateParameter("EMIS.XYL", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.4e-11));
  state.SetCustomRateParameter("EMIS.ISOP", std::vector<double>(solver.parameters_.number_of_grid_cells_, 6.03e-12));
  state.SetCustomRateParameter("PHOTO.NO2", std::vector<double>(solver.parameters_.number_of_grid_cells_, 0.00477));
  state.SetCustomRateParameter("PHOTO.O3->O1D", std::vector<double>(solver.parameters_.number_of_grid_cells_, 2.26e-06));
  state.SetCustomRateParameter("PHOTO.O3->O3P", std::vector<double>(solver.parameters_.number_of_grid_cells_, 0.000253));
  state.SetCustomRateParameter("PHOTO.NO3->NO2", std::vector<double>(solver.parameters_.number_of_grid_cells_, 0.117));
  state.SetCustomRateParameter("PHOTO.NO3->NO", std::vector<double>(solver.parameters_.number_of_grid_cells_, 0.0144));
  state.SetCustomRateParameter("PHOTO.HONO", std::vector<double>(solver.parameters_.number_of_grid_cells_, 0.000918));
  state.SetCustomRateParameter("PHOTO.H2O2", std::vector<double>(solver.parameters_.number_of_grid_cells_, 2.59e-06));
  state.SetCustomRateParameter("PHOTO.PNA", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.89e-06));
  state.SetCustomRateParameter("PHOTO.HNO3", std::vector<double>(solver.parameters_.number_of_grid_cells_, 8.61e-08));
  state.SetCustomRateParameter("PHOTO.NTR", std::vector<double>(solver.parameters_.number_of_grid_cells_, 4.77e-07));
  state.SetCustomRateParameter("PHOTO.ROOH", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.81e-06));
  state.SetCustomRateParameter("PHOTO.MEPX", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.81e-06));
  state.SetCustomRateParameter("PHOTO.FORM->HO2", std::vector<double>(solver.parameters_.number_of_grid_cells_, 7.93e-06));
  state.SetCustomRateParameter("PHOTO.FORM->CO", std::vector<double>(solver.parameters_.number_of_grid_cells_, 2.2e-05));
  state.SetCustomRateParameter("PHOTO.ALD2", std::vector<double>(solver.parameters_.number_of_grid_cells_, 2.2e-06));
  state.SetCustomRateParameter("PHOTO.PACD", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.81e-06));
  state.SetCustomRateParameter("PHOTO.ALDX", std::vector<double>(solver.parameters_.number_of_grid_cells_, 2.2e-06));
  state.SetCustomRateParameter("PHOTO.OPEN", std::vector<double>(solver.parameters_.number_of_grid_cells_, 0.000645));
  state.SetCustomRateParameter("PHOTO.MGLY", std::vector<double>(solver.parameters_.number_of_grid_cells_, 7.64e-05));
  state.SetCustomRateParameter("PHOTO.ISPD", std::vector<double>(solver.parameters_.number_of_grid_cells_, 1.98e-09));

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
      auto result = solver.template Solve<true>(time_step - elapsed_solve_time, state);
      auto end = std::chrono::high_resolution_clock::now();
      total_solve_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
      elapsed_solve_time = result.final_time_;
      state.variables_ = result.result_;

      total_stats.function_calls += result.stats_.function_calls;
      total_stats.jacobian_updates += result.stats_.jacobian_updates;
      total_stats.number_of_steps += result.stats_.number_of_steps;
      total_stats.accepted += result.stats_.accepted;
      total_stats.rejected += result.stats_.rejected;
      total_stats.decompositions += result.stats_.decompositions;
      total_stats.solves += result.stats_.solves;
      total_stats.singular += result.stats_.singular;
      total_stats.total_update_state_time += result.stats_.total_update_state_time;
      total_stats.total_forcing_time += result.stats_.total_forcing_time;
      total_stats.total_jacobian_time += result.stats_.total_jacobian_time;
      total_stats.total_linear_factor_time += result.stats_.total_linear_factor_time;
      total_stats.total_linear_solve_time += result.stats_.total_linear_solve_time;
    }
  }

  return std::make_tuple(state, total_stats, total_solve_time);
}

int main(const int argc, const char* argv[])
{
  SolverConfig solverConfig;

  std::string config_path = "./configs/carbon_bond_5";
  ConfigParseStatus status = solverConfig.ReadAndParse(config_path);
  if (status != micm::ConfigParseStatus::Success)
  {
    throw "Parsing failed";
  }

  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  auto chemical_system = solver_params.system_;
  auto reactions = solver_params.processes_;

  auto solver_parameters = RosenbrockSolverParameters::three_stage_rosenbrock_parameters(n_grid_cells);

  RosenbrockSolver<GroupVectorMatrix, GroupSparseVectorMatrix> solver{ chemical_system, reactions, solver_parameters };

  auto jit{ micm::JitCompiler::create() };

  auto start = std::chrono::high_resolution_clock::now();
  JitRosenbrockSolver<
      GroupVectorMatrix,
      GroupSparseVectorMatrix,
      JitLinearSolver<n_grid_cells, GroupSparseVectorMatrix>,
      JitProcessSet<n_grid_cells>>
      jit_solver(jit.get(), chemical_system, reactions, solver_parameters);
  auto end = std::chrono::high_resolution_clock::now();
  auto jit_compile_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

  std::cout << "Jit compile time: " << jit_compile_time.count() << " nanoseconds" << std::endl;

  auto result_tuple = run_solver(solver);
  auto jit_result_tuple = run_solver(jit_solver);

  std::cout << "Standard solve time: " << std::get<2>(result_tuple).count() << " nanoseconds" << std::endl;
  std::cout << "JIT solve time: " << std::get<2>(jit_result_tuple).count() << " nanoseconds" << std::endl;

  auto result_stats = std::get<1>(result_tuple);
  std::cout << "Standard solve stats: " << std::endl;
  std::cout << "\taccepted: " << result_stats.accepted << std::endl;
  std::cout << "\tfunction_calls: " << result_stats.function_calls << std::endl;
  std::cout << "\tjacobian_updates: " << result_stats.jacobian_updates << std::endl;
  std::cout << "\tnumber_of_steps: " << result_stats.number_of_steps << std::endl;
  std::cout << "\taccepted: " << result_stats.accepted << std::endl;
  std::cout << "\trejected: " << result_stats.rejected << std::endl;
  std::cout << "\tdecompositions: " << result_stats.decompositions << std::endl;
  std::cout << "\tsolves: " << result_stats.solves << std::endl;
  std::cout << "\tsingular: " << result_stats.singular << std::endl;
  std::cout << "\ttotal_update_state_time: " << result_stats.total_update_state_time.count() << " nanoseconds" << std::endl;
  std::cout << "\ttotal_forcing_time: " << result_stats.total_forcing_time.count() << " nanoseconds" << std::endl;
  std::cout << "\ttotal_jacobian_time: " << result_stats.total_jacobian_time.count() << " nanoseconds" << std::endl;
  std::cout << "\ttotal_linear_factor_time: " << result_stats.total_linear_factor_time.count() << " nanoseconds" << std::endl;
  std::cout << "\ttotal_linear_solve_time: " << result_stats.total_linear_solve_time.count() << " nanoseconds" << std::endl << std::endl;

  auto jit_result_stats = std::get<1>(jit_result_tuple);
  std::cout << "JIT solve stats: " << std::endl;
  std::cout << "\taccepted: " << jit_result_stats.accepted << std::endl;
  std::cout << "\tfunction_calls: " << jit_result_stats.function_calls << std::endl;
  std::cout << "\tjacobian_updates: " << jit_result_stats.jacobian_updates << std::endl;
  std::cout << "\tnumber_of_steps: " << jit_result_stats.number_of_steps << std::endl;
  std::cout << "\taccepted: " << jit_result_stats.accepted << std::endl;
  std::cout << "\trejected: " << jit_result_stats.rejected << std::endl;
  std::cout << "\tdecompositions: " << jit_result_stats.decompositions << std::endl;
  std::cout << "\tsolves: " << jit_result_stats.solves << std::endl;
  std::cout << "\tsingular: " << jit_result_stats.singular << std::endl;
  std::cout << "\ttotal_update_state_time: " << jit_result_stats.total_update_state_time.count() << " nanoseconds" << std::endl;
  std::cout << "\ttotal_forcing_time: " << jit_result_stats.total_forcing_time.count() << " nanoseconds" << std::endl;
  std::cout << "\ttotal_jacobian_time: " << jit_result_stats.total_jacobian_time.count() << " nanoseconds" << std::endl;
  std::cout << "\ttotal_linear_factor_time: " << jit_result_stats.total_linear_factor_time.count() << " nanoseconds"
            << std::endl;
  std::cout << "\ttotal_linear_solve_time: " << jit_result_stats.total_linear_solve_time.count() << " nanoseconds"
            << std::endl;

  auto result = std::get<0>(result_tuple);
  auto jit_result = std::get<0>(jit_result_tuple);

  for (auto& species : result.variable_names_)
  {
    for (int i = 0; i < n_grid_cells; ++i)
    {
      double a = result.variables_[i][result.variable_map_[species]];
      double b = jit_result.variables_[i][jit_result.variable_map_[species]];
      if ( std::abs(a - b) > 1.0e-5 * (std::abs(a) + std::abs(b)) / 2.0 + 1.0e-30 )
      {
        std::cout << species << " does not match final concentration" << std::endl;
      }
    }
  }
  return 0;
}