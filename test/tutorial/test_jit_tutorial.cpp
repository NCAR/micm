#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>

// Each rate constant is in its own header file
#include <micm/configure/solver_config.hpp>
#include <micm/jit/jit_compiler.hpp>
#include <micm/solver/jit_rosenbrock.hpp>
#include <micm/solver/rosenbrock.hpp>

// Use our namespace so that this example is easier to read
using namespace micm;

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;

template<typename T>
auto run_solver(T& solver)
{
  typename T::SolverStats total_stats;
  State state = solver.GetState();

  state.variables_ = 1;

  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa

  state.SetCustomRateParameter("EMIS.NO", 1.44e-10);
  state.SetCustomRateParameter("EMIS.NO2", 7.56e-12);
  state.SetCustomRateParameter("EMIS.CO", 1.9600000000000003e-09);
  state.SetCustomRateParameter("EMIS.SO2", 1.06e-09);
  state.SetCustomRateParameter("EMIS.FORM", 1.02e-11);
  state.SetCustomRateParameter("EMIS.MEOH", 5.920000000000001e-13);
  state.SetCustomRateParameter("EMIS.ALD2", 4.25e-12);
  state.SetCustomRateParameter("EMIS.PAR", 4.27e-10);
  state.SetCustomRateParameter("EMIS.ETH", 4.62e-11);
  state.SetCustomRateParameter("EMIS.OLE", 1.49e-11);
  state.SetCustomRateParameter("EMIS.IOLE", 1.49e-11);
  state.SetCustomRateParameter("EMIS.TOL", 1.53e-11);
  state.SetCustomRateParameter("EMIS.XYL", 1.4e-11);
  state.SetCustomRateParameter("EMIS.ISOP", 6.03e-12);
  state.SetCustomRateParameter("PHOTO.NO2", 0.00477);
  state.SetCustomRateParameter("PHOTO.O3->O1D", 2.26e-06);
  state.SetCustomRateParameter("PHOTO.O3->O3P", 0.00025299999999999997);
  state.SetCustomRateParameter("PHOTO.NO3->NO2", 0.11699999999999999);
  state.SetCustomRateParameter("PHOTO.NO3->NO", 0.0144);
  state.SetCustomRateParameter("PHOTO.HONO", 0.000918);
  state.SetCustomRateParameter("PHOTO.H2O2", 2.59e-06);
  state.SetCustomRateParameter("PHOTO.PNA", 1.89e-06);
  state.SetCustomRateParameter("PHOTO.HNO3", 8.61e-08);
  state.SetCustomRateParameter("PHOTO.NTR", 4.77e-07);
  state.SetCustomRateParameter("PHOTO.ROOH", 1.81e-06);
  state.SetCustomRateParameter("PHOTO.MEPX", 1.81e-06);
  state.SetCustomRateParameter("PHOTO.FORM->HO2", 7.93e-06);
  state.SetCustomRateParameter("PHOTO.FORM->CO", 2.2e-05);
  state.SetCustomRateParameter("PHOTO.ALD2", 2.2e-06);
  state.SetCustomRateParameter("PHOTO.PACD", 1.81e-06);
  state.SetCustomRateParameter("PHOTO.ALDX", 2.2e-06);
  state.SetCustomRateParameter("PHOTO.OPEN", 0.0006450000000000001);
  state.SetCustomRateParameter("PHOTO.MGLY", 7.64e-05);
  state.SetCustomRateParameter("PHOTO.ISPD", 1.98e-09);

  // choose a timestep and print the initial state
  double time_step = 500;  // s

  // solve for ten iterations
  for (int i = 0; i < 10; ++i)
  {
    double elapsed_solve_time = 0;

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.template Solve<true>(time_step - elapsed_solve_time, state);
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
      total_stats.total_forcing_time += result.stats_.total_forcing_time;
      total_stats.total_jacobian_time += result.stats_.total_jacobian_time;
      total_stats.total_linear_factor_time += result.stats_.total_linear_factor_time;
      total_stats.total_linear_solve_time += result.stats_.total_linear_solve_time;
    }
  }

  return std::make_pair(state, total_stats);
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

  RosenbrockSolver<> solver{ chemical_system, reactions, RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };

  auto jit{ micm::JitCompiler::create() };

  auto start = std::chrono::high_resolution_clock::now();
  JitRosenbrockSolver<
    Group1VectorMatrix, 
    Group1SparseVectorMatrix, 
    JitLinearSolver<1, Group1SparseVectorMatrix>
    > jit_solver(
      jit.get(), 
      chemical_system, 
      reactions, 
      RosenbrockSolverParameters::three_stage_rosenbrock_parameters() 
    );
  auto end = std::chrono::high_resolution_clock::now();
  auto jit_compile_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

  std::cout << "Jit compile time: " << jit_compile_time.count() << " nanoseconds" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  auto result_pair = run_solver(solver);
  end = std::chrono::high_resolution_clock::now();
  auto result_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

  start = std::chrono::high_resolution_clock::now();
  auto jit_result_pair = run_solver(jit_solver);
  end = std::chrono::high_resolution_clock::now();
  auto jit_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

  std::cout << "Result time: " << result_time.count() << " nanoseconds" << std::endl;
  std::cout << "JIT result time: " << jit_time.count() << " nanoseconds" << std::endl;

  auto result_stats = result_pair.second;
  std::cout << "Result stats: " << std::endl;
  std::cout << "accepted: " << result_stats.accepted << std::endl;
  std::cout << "function_calls: " << result_stats.function_calls << std::endl;
  std::cout << "jacobian_updates: " << result_stats.jacobian_updates << std::endl;
  std::cout << "number_of_steps: " << result_stats.number_of_steps << std::endl;
  std::cout << "accepted: " << result_stats.accepted << std::endl;
  std::cout << "rejected: " << result_stats.rejected << std::endl;
  std::cout << "decompositions: " << result_stats.decompositions << std::endl;
  std::cout << "solves: " << result_stats.solves << std::endl;
  std::cout << "singular: " << result_stats.singular << std::endl;
  std::cout << "total_forcing_time: " <<       result_stats.total_forcing_time.count() << " nanoseconds" << std::endl;
  std::cout << "total_jacobian_time: " <<      result_stats.total_jacobian_time.count() << " nanoseconds" << std::endl;
  std::cout << "total_linear_factor_time: " << result_stats.total_linear_factor_time.count() << " nanoseconds" << std::endl;
  std::cout << "total_linear_solve_time: " <<  result_stats.total_linear_solve_time.count() << " nanoseconds" << std::endl;

  auto jit_result_stats = jit_result_pair.second;
  std::cout << "JIT result stats: " << std::endl;
  std::cout << "accepted: " << jit_result_stats.accepted << std::endl;
  std::cout << "function_calls: " << jit_result_stats.function_calls << std::endl;
  std::cout << "jacobian_updates: " << jit_result_stats.jacobian_updates << std::endl;
  std::cout << "number_of_steps: " << jit_result_stats.number_of_steps << std::endl;
  std::cout << "accepted: " << jit_result_stats.accepted << std::endl;
  std::cout << "rejected: " << jit_result_stats.rejected << std::endl;
  std::cout << "decompositions: " << jit_result_stats.decompositions << std::endl;
  std::cout << "solves: " << jit_result_stats.solves << std::endl;
  std::cout << "singular: " << jit_result_stats.singular << std::endl;
  std::cout << "total_forcing_time: " <<       jit_result_stats.total_forcing_time.count() << " nanoseconds" << std::endl;
  std::cout << "total_jacobian_time: " <<      jit_result_stats.total_jacobian_time.count() << " nanoseconds" << std::endl;
  std::cout << "total_linear_factor_time: " << jit_result_stats.total_linear_factor_time.count() << " nanoseconds" << std::endl;
  std::cout << "total_linear_solve_time: " <<  jit_result_stats.total_linear_solve_time.count() << " nanoseconds" << std::endl;

  auto result = result_pair.first;
  auto jit_result = result_pair.first;

  for(auto& species : result.variable_names_) {
    if (result.variables_[0][result.variable_map_[species]] != jit_result.variables_[0][jit_result.variable_map_[species]]) {
      std::cout << species << " do not match final concentration" << std::endl;
    }
  }
  return 0;
}
