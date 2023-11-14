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

  RosenbrockSolver<> solver{ System(SystemParameters{ .gas_phase_ = gas_phase }),
                             std::vector<Process>{ r1, r2, r3 },
                             params };

  auto state = solver.GetState();

  // mol m-3
  state.SetConcentration(a, std::vector<double>(3, 1));
  state.SetConcentration(b, std::vector<double>(3, 2));
  state.SetConcentration(c, std::vector<double>(3, 3));

  for (auto& elem : state.variables_.AsVector())
  {
    std::cout << elem << std::endl;
  }

  std::cout << std::endl;

  RosenbrockSolver<Group3Matrix, Group3SparseVectorMatrix> vectorized_solver{
    System(SystemParameters{ .gas_phase_ = gas_phase }), std::vector<Process>{ r1, r2, r3 }, params
  };

  auto vectorized_state = solver.GetState();

  // mol m-3
  vectorized_state.SetConcentration(a, std::vector<double>(3, 1));
  vectorized_state.SetConcentration(b, std::vector<double>(3, 2));
  vectorized_state.SetConcentration(c, std::vector<double>(3, 3));

  for (auto& elem : vectorized_state.variables_.AsVector())
  {
    std::cout << elem << std::endl;
  }

  // double k1 = 0.04;
  // double k2 = 3e7;
  // double k3 = 1e4;
  // state.SetCustomRateParameter("r1", std::vector<double>(3, k1));
  // state.SetCustomRateParameter("r2", std::vector<double>(3, k2));
  // state.SetCustomRateParameter("r3", std::vector<double>(3, k3));

  // double temperature = 272.5;  // [K]
  // double pressure = 101253.3;  // [Pa]
  // double air_density = 1e6;    // [mol m-3]

  // for (size_t cell = 0; cell < solver.parameters_.number_of_grid_cells_; ++cell)
  // {
  //   state.conditions_[cell].temperature_ = temperature;
  //   state.conditions_[cell].pressure_ = pressure;
  //   state.conditions_[cell].air_density_ = air_density;
  // }

  // // choose a timestep and print the initial state
  // double time_step = 200;  // s

  // auto result = solver.Solve(time_step, state);

  // for (int i = 0; i < 10; ++i)
  // {
  //   double elapsed_solve_time = 0;
  //   while (elapsed_solve_time < time_step)
  //   {
  //     auto result = solver.Solve(time_step - elapsed_solve_time, state);
  //     elapsed_solve_time = result.final_time_;
  //     state.variables_[0] = result.result_.AsVector();
  //   }

  // }
}