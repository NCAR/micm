#include <chrono>
#include <iomanip>
#include <iostream>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

// Use our namespace so that this example is easier to read
using namespace micm;

template<typename T>
using Group3Matrix = micm::VectorMatrix<T, 3>;

template<typename T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;

void solve(auto& solver, auto& state)
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

  for (size_t cell = 0; cell < solver.parameters_.number_of_grid_cells_; ++cell)
  {
    state.conditions_[cell].temperature_ = temperature;
    state.conditions_[cell].pressure_ = pressure;
    state.conditions_[cell].air_density_ = air_density;
  }

  double time_step = 100;  // s

  auto result = solver.Solve(time_step, state);

  for (int i = 0; i < 10; ++i)
  {
    double elapsed_solve_time = 0;
    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
      state.variables_ = result.result_;
    }
  }
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

  auto params = RosenbrockSolverParameters::three_stage_rosenbrock_parameters(3, false);
  auto system = System(SystemParameters{ .gas_phase_ = gas_phase });
  auto reactions = std::vector<Process>{ r1, r2, r3 };

  RosenbrockSolver<> solver{ system, reactions, params };

  auto state = solver.GetState();

  state.SetConcentration(a, { 1.1, 2.1, 3.1 });
  state.SetConcentration(b, { 1.2, 2.2, 3.2 });
  state.SetConcentration(c, { 1.3, 2.3, 3.3 });

  for (auto& elem : state.variables_.AsVector())
  {
    std::cout << elem << std::endl;
  }

  std::cout << std::endl;

  RosenbrockSolver<Group3Matrix, Group3SparseVectorMatrix> vectorized_solver{ system, reactions, params };

  auto vectorized_state = vectorized_solver.GetState();

  vectorized_state.SetConcentration(a, { 1.1, 2.1, 3.1 });
  vectorized_state.SetConcentration(b, { 1.2, 2.2, 3.2 });
  vectorized_state.SetConcentration(c, { 1.3, 2.3, 3.3 });

  for (auto& elem : vectorized_state.variables_.AsVector())
  {
    std::cout << elem << std::endl;
  }

  std::cout << std::endl;

  solve(solver, state);
  solve(vectorized_solver, vectorized_state);

  for (size_t cell = 0; cell < params.number_of_grid_cells_; ++cell)
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

  micm::SparseMatrix<std::string> sparse_matrix{ micm::SparseMatrix<std::string>::create(4)
                                                     .with_element(0, 1)
                                                     .with_element(2, 1)
                                                     .with_element(2, 3)
                                                     .with_element(3, 2)
                                                     .number_of_blocks(3) };

  micm::SparseMatrix<std::string, micm::SparseMatrixVectorOrdering<3>> sparse_vector_matrix{
    micm::SparseMatrix<std::string, micm::SparseMatrixVectorOrdering<3>>::create(4)
        .with_element(0, 1)
        .with_element(2, 1)
        .with_element(2, 3)
        .with_element(3, 2)
        .number_of_blocks(3)
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