#include <iomanip>
#include <iostream>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

using namespace micm;

template<class T>
using SparseMatrixPolicy = SparseMatrix<T>;

int main(const int argc, const char *argv[])
{
  auto foo = Species{ "Foo" };
  auto bar = Species{ "Bar" };
  auto baz = Species{ "Baz" };

  Phase gas_phase{ std::vector<Species>{ foo, bar, baz } };

  System chemical_system{ SystemParameters{ .gas_phase_ = gas_phase } };

  Process r1 = Process::create()
                   .reactants({ foo })
                   .products({ Yield(bar, 0.8), Yield(baz, 0.2) })
                   .rate_constant(ArrheniusRateConstant({ .A_ = 1.0e-3 }))
                   .phase(gas_phase);

  Process r2 = Process::create()
                   .reactants({ foo, bar })
                   .products({ Yield(baz, 1) })
                   .rate_constant(ArrheniusRateConstant({ .A_ = 1.0e-5, .C_ = 110.0 }))
                   .phase(gas_phase);

  std::vector<Process> reactions{ r1, r2 };

  RosenbrockSolver<Matrix, SparseMatrixPolicy> solver{ chemical_system, reactions, RosenbrockSolverParameters{} };

  State state = solver.GetState();

  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa
  state.SetConcentration(foo, 20.0);           // mol m-3

  std::cout << "foo, bar, baz" << std::endl;
  for (int i = 0; i < 10; ++i)
  {
    auto result = solver.Solve(500.0, state);
    state.variables_ = result.result_;
    std::cout << std::fixed << std::setprecision(6) << state.variables_[0][state.variable_map_["Foo"]] << ", "
              << state.variables_[0][state.variable_map_["Bar"]] << ", " << state.variables_[0][state.variable_map_["Baz"]]
              << std::endl;
  }

  return 0;
}