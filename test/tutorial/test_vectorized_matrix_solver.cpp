#include <micm/CPU.hpp>

#include <chrono>
#include <iomanip>
#include <iostream>

// Use our namespace so that this example is easier to read
using namespace micm;

void solve(auto& solver, auto& state, std::size_t number_of_grid_cells)
{
  double k1 = 0.04;
  double k2 = 3e7;
  double k3 = 1e4;
  state.SetCustomRateParameter("r1", std::vector<double>(3, k1));
  state.SetCustomRateParameter("r2", std::vector<double>(3, k2));
  state.SetCustomRateParameter("r3", std::vector<double>(3, k3));

  double temperature = 272.5;  // [K]
  double pressure = 101253.3;  // [Pa]
  double air_density = 1e6;    // [mol m-3]

  for (size_t cell = 0; cell < number_of_grid_cells; ++cell)
  {
    state.conditions_[cell].temperature_ = temperature;
    state.conditions_[cell].pressure_ = pressure;
    state.conditions_[cell].air_density_ = air_density;
  }

  double time_step = 100;  // s

  SolverState solver_state = SolverState::Converged;
  for (int i = 0; i < 10 && solver_state == SolverState::Converged; ++i)
  {
    double elapsed_solve_time = 0;
    solver.CalculateRateConstants(state);
    while (elapsed_solve_time < time_step && solver_state != SolverState::Converged)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time += result.stats_.final_time_;;
      solver_state = result.state_;
    }
  }
  if (solver_state != SolverState::Converged)
  {
    throw "Solver did not converge";
  }
}

int main()
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species("C");

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ a, b, c } };

  Process r1 = ChemicalReactionBuilder()
                   .SetReactants({ a })
                   .SetProducts({ Yield(b, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "r1" }))
                   .SetPhase(gas_phase)
                   .Build();

  Process r2 = ChemicalReactionBuilder()
                   .SetReactants({ b, b })
                   .SetProducts({ Yield(b, 1), Yield(c, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "r2" }))
                   .SetPhase(gas_phase)
                   .Build();

  Process r3 = ChemicalReactionBuilder()
                   .SetReactants({ b, c })
                   .SetProducts({ Yield(a, 1), Yield(c, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "r3" }))
                   .SetPhase(gas_phase)
                   .Build();

  auto params = RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  auto system = System(SystemParameters{ .gas_phase_ = gas_phase });
  auto reactions = std::vector<Process>{ r1, r2, r3 };
  const std::size_t number_of_grid_cells = 3;

  auto solver = CpuSolverBuilder<micm::RosenbrockSolverParameters>(params).SetSystem(system).SetReactions(reactions).Build();

  auto state = solver.GetState(number_of_grid_cells);

  state.SetConcentration(a, { 1.1, 2.1, 3.1 });
  state.SetConcentration(b, { 1.2, 2.2, 3.2 });
  state.SetConcentration(c, { 1.3, 2.3, 3.3 });

  for (auto& elem : state.variables_.AsVector())
  {
    std::cout << elem << std::endl;
  }

  std::cout << std::endl;

  auto vectorized_solver = CpuSolverBuilder<
                               RosenbrockSolverParameters,
                               VectorMatrix<double, 3>,
                               SparseMatrix<double, SparseMatrixVectorOrdering<3>>>(params)
                               .SetSystem(system)
                               .SetReactions(reactions)
                               .Build();

  auto vectorized_state = vectorized_solver.GetState(number_of_grid_cells);

  vectorized_state.SetConcentration(a, { 1.1, 2.1, 3.1 });
  vectorized_state.SetConcentration(b, { 1.2, 2.2, 3.2 });
  vectorized_state.SetConcentration(c, { 1.3, 2.3, 3.3 });

  for (auto& elem : vectorized_state.variables_.AsVector())
  {
    std::cout << elem << std::endl;
  }

  std::cout << std::endl;

  solve(solver, state, number_of_grid_cells);
  solve(vectorized_solver, vectorized_state, number_of_grid_cells);

  for (size_t cell = 0; cell < number_of_grid_cells; ++cell)
  {
    std::cout << "Cell " << cell << std::endl;
    std::cout << std::setw(10) << "Species" << std::setw(20) << "Regular" << std::setw(20) << "Vectorized" << std::endl;

    for (auto& species : system.UniqueNames())
    {
      auto regular_idx = state.variable_map_[species];
      auto vectorized_idx = vectorized_state.variable_map_[species];

      std::cout << std::setw(10) << species << std::setw(20) << state.variables_[cell][regular_idx] << std::setw(20)
                << vectorized_state.variables_[cell][vectorized_idx] << std::endl;
    }

    std::cout << std::endl;
  }

  micm::SparseMatrix<std::string> sparse_matrix{ micm::SparseMatrix<std::string>::Create(4)
                                                     .WithElement(0, 1)
                                                     .WithElement(2, 1)
                                                     .WithElement(2, 3)
                                                     .WithElement(3, 2)
                                                     .SetNumberOfBlocks(3) };

  micm::SparseMatrix<std::string, micm::SparseMatrixVectorOrdering<3>> sparse_vector_matrix{
    micm::SparseMatrix<std::string, micm::SparseMatrixVectorOrdering<3>>::Create(4)
        .WithElement(0, 1)
        .WithElement(2, 1)
        .WithElement(2, 3)
        .WithElement(3, 2)
        .SetNumberOfBlocks(3)
  };

  sparse_matrix[0][0][1] = sparse_vector_matrix[0][0][1] = "0.0.1";
  sparse_matrix[0][2][1] = sparse_vector_matrix[0][2][1] = "0.2.1";
  sparse_matrix[0][2][3] = sparse_vector_matrix[0][2][3] = "0.2.3";
  sparse_matrix[0][3][2] = sparse_vector_matrix[0][3][2] = "0.3.2";
  sparse_matrix[1][0][1] = sparse_vector_matrix[1][0][1] = "1.0.1";
  sparse_matrix[1][2][1] = sparse_vector_matrix[1][2][1] = "1.2.1";
  sparse_matrix[1][2][3] = sparse_vector_matrix[1][2][3] = "1.2.3";
  sparse_matrix[1][3][2] = sparse_vector_matrix[1][3][2] = "1.3.2";
  sparse_matrix[2][0][1] = sparse_vector_matrix[2][0][1] = "2.0.1";
  sparse_matrix[2][2][1] = sparse_vector_matrix[2][2][1] = "2.2.1";
  sparse_matrix[2][2][3] = sparse_vector_matrix[2][2][3] = "2.2.3";
  sparse_matrix[2][3][2] = sparse_vector_matrix[2][3][2] = "2.3.2";

  std::cout << std::endl << "Sparse matrix standard ordering elements" << std::endl;
  for (auto& elem : sparse_matrix.AsVector())
  {
    std::cout << elem << std::endl;
  }
  std::cout << std::endl << "Sparse matrix vector ordering elements" << std::endl;
  for (auto& elem : sparse_vector_matrix.AsVector())
  {
    std::cout << elem << std::endl;
  }
}
