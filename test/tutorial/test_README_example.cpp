#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>

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

  auto solver = micm::CpuSolverBuilder(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(chemical_system)
                    .SetReactions(reactions)
                    .Build();

  State state = solver.GetState();

  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa
  state.SetConcentration(foo, 20.0);           // mol m-3

  state.PrintHeader();
  for (int i = 0; i < 10; ++i)
  {
    auto result = solver.Solve(500.0, state);
    state.PrintState(i * 500);
  }

  return 0;
}