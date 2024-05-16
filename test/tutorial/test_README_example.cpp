#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

#include <iomanip>
#include <iostream>

using namespace micm;

int main(const int argc, const char *argv[])
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

#if 0
  auto rosenbrock_solver = Solver::Create()
                    .SetSystem(chemical_system)
                    .SetReactions(reactions)
                    .SetNumberOfGridCells(19)
                    .Rosenbrock<MatrixPolicy, SparseMatrixPolicy>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters());

  auto backward_euler_solver = Solver::Create()
                    .SetSystem(chemical_system)
                    .SetReactions(reactions)
                    .SetNumberOfGridCells(19)
                    .BackwardEuler(BackwardEulerParameters::CamChemParameters());

  auto cuda_solver = Solver::Create()
                    .SetSystem(chemical_system)
                    .SetReactions(reactions)
                    .SetNumberOfGridCells(19)
                    .Rosenbrock<CudaDenseMatrixPolicy, CudaSparseMatrixPolicy>(RosenbrockSolverParameters::ThreeStageRosenbrockParameter());
  

  auto rosenbrock_state = rosenbrock_solver.GetState();
  auto backward_euler_state = backward_euler_solver.GetState();
  auto cuda_state = cuda_solver.GetState();

  rosenbrock_state.conditions_[0].temperature_ = 287.45;  // K
  rosenbrock_state.conditions_[0].pressure_ = 101319.9;   // Pa
  rosenbrock_state.SetConcentration(foo, 20.0);           // mol m-3

  backward_euler_state.conditions_[0].temperature_ = 287.45;  // K
  backward_euler_state.conditions_[0].pressure_ = 101319.9;   // Pa
  backward_euler_state.SetConcentration(foo, 20.0);           // mol m-3

  cuda_state.conditions_[0].temperature_ = 287.45;  // K
  cuda_state.conditions_[0].pressure_ = 101319.9;   // Pa
  cuda_state.SetConcentration(foo, 20.0);           // mol m-3

  cuda_state.PrintHeader();
  for (int i = 0; i < 10; ++i)
  {
    cuda_state.CalculateRateConstants(); // calculate rate constants on the host
    cuda_state.SyncToDevice();
    auto result = cuda_solver.Solve(500.0, cuda_state); // solve on the device
    cuda_state.variables_ = result.result_; // device to device copy with overloaded operator=
    cuda_state.SyncToHost();
    cuda_state.PrintState(i * 500);
  }
#endif
  return 0;
}