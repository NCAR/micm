// Each rate constant is in its own header file
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/solver_builder.hpp>

#include <iomanip>
#include <iostream>
#include <map>

// Use our namespace so that this example is easier to read
using namespace micm;

int main(const int argc, const char* argv[])
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species(
      "C",
      std::map<std::string, double>{ { "molecular weight [kg mol-1]", 0.025 },
                                     { "diffusion coefficient [m2 s-1]", 2.3e2 } });
  auto d = Species("D");
  auto e = Species("E");
  auto f = Species("F");
  auto g = Species("G");

  Phase gas_phase{ std::vector<Species>{ a, b, c, d, e, f, g } };

  Process r1 = Process::Create()
                   .SetReactants({ a })
                   .SetProducts({ Yields(b, 1) })
                   .SetRateConstant(ArrheniusRateConstant({ .A_ = 2.15e-1, .B_ = 0, .C_ = 110 }))
                   .SetPhase(gas_phase);

  // a branched reaction has two output pathways
  // this is represnted internal to micm as two different reactions
  auto branched_params = BranchedRateConstantParameters{ .X_ = 1.2, .Y_ = 204.3, .a0_ = 1.0e-3, .n_ = 2 };
  branched_params.branch_ = BranchedRateConstantParameters::Branch::Alkoxy;

  Process r2 = Process::Create()
                   .SetReactants({ b })
                   .SetProducts({ Yields(c, 1) })
                   .SetRateConstant(BranchedRateConstant(branched_params))
                   .SetPhase(gas_phase);

  branched_params.branch_ = BranchedRateConstantParameters::Branch::Nitrate;
  Process r3 = Process::Create()
                   .SetReactants({ b })
                   .SetProducts({ Yields(d, 1) })
                   .SetRateConstant(BranchedRateConstant(branched_params))
                   .SetPhase(gas_phase);

  // A surface rate constant also needs to know the effective radius and particle number concentration
  // we will set those later
  Process r4 = Process::Create()
                   .SetReactants({ c })
                   .SetProducts({ Yields(e, 1) })
                   .SetRateConstant(SurfaceRateConstant({ .label_ = "C", .species_ = c, .reaction_probability_ = 0.90 }))
                   .SetPhase(gas_phase);

  Process r5 = Process::Create()
                   .SetReactants({ d })
                   .SetProducts({ Yields(f, 2) })
                   .SetRateConstant(TernaryChemicalActivationRateConstant({ .k0_A_ = 1.2,
                                                                          .k0_B_ = 2.3,
                                                                          .k0_C_ = 302.3,
                                                                          .kinf_A_ = 2.6,
                                                                          .kinf_B_ = -3.1,
                                                                          .kinf_C_ = 402.1,
                                                                          .Fc_ = 0.9,
                                                                          .N_ = 1.2 }))
                   .SetPhase(gas_phase);

  // to have a stoichiemetric coefficient of more than one for reactants,
  // list the reactant that many times
  Process r6 = Process::Create()
                   .SetReactants({ e, e })
                   .SetProducts({ Yields(g, 1) })
                   .SetRateConstant(TroeRateConstant({ .k0_A_ = 1.2e4,
                                                     .k0_B_ = 167.0,
                                                     .k0_C_ = 3.0,
                                                     .kinf_A_ = 136.0,
                                                     .kinf_B_ = 5.0,
                                                     .kinf_C_ = 24.0,
                                                     .Fc_ = 0.9,
                                                     .N_ = 0.8 }))
                   .SetPhase(gas_phase);

  Process r7 = Process::Create()
                   .SetReactants({ f })
                   .SetProducts({ Yields(g, 1) })
                   .SetRateConstant(TunnelingRateConstant({ .A_ = 1.2, .B_ = 2.3, .C_ = 302.3 }))
                   .SetPhase(gas_phase);

  Process r8 = Process::Create()
                   .SetReactants({ c })
                   .SetProducts({ Yields(g, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "my photolysis rate" }))
                   .SetPhase(gas_phase);

  Process r9 = Process::Create()
                   .SetProducts({ Yields(a, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "my emission rate" }))
                   .SetPhase(gas_phase);

  Process r10 = Process::Create()
                    .SetReactants({ b })
                    .SetRateConstant(UserDefinedRateConstant({ .label_ = "my loss rate" }))
                    .SetPhase(gas_phase);

  auto chemical_system = System(micm::SystemParameters{ .gas_phase_ = gas_phase });
  auto reactions = std::vector<micm::Process>{ r1, r2, r3, r4, r5, r6, r7, r8, r9, r10 };

  auto solver = micm::CpuSolverBuilder_DoolittleLU<micm::RosenbrockSolverParameters>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(chemical_system)
                    .SetReactions(reactions)
                    .Build();
  State state = solver.GetState();

  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa

  state.SetConcentration(a, 1.0);  // mol m-3
  state.SetConcentration(b, 0.0);  // mol m-3
  state.SetConcentration(c, 0.0);  // mol m-3
  state.SetConcentration(d, 0.0);  // mol m-3
  state.SetConcentration(e, 0.0);  // mol m-3
  state.SetConcentration(f, 0.0);  // mol m-3
  state.SetConcentration(g, 0.0);  // mol m-3

  state.SetCustomRateParameter("C.effective radius [m]", 1e-7);
  state.SetCustomRateParameter("C.particle number concentration [# m-3]", 2.5e6);

  // choose and timestep a print the initial state
  double time_step = 500;  // s

  state.PrintHeader();
  state.PrintState(0);

  double photo_rate = 1e-10;
  double emission_rate = 1e-20;
  double loss = emission_rate * 1e-3;
  // these rates are constant through the simulation
  state.SetCustomRateParameter("my emission rate", emission_rate);
  state.SetCustomRateParameter("my loss rate", loss);

  // solve for ten iterations
  for (int i = 0; i < 10; ++i)
  {
    // Depending on how stiff the system is
    // the solver integration step may not be able to solve for the full time step
    // so we need to track how much time the solver was able to integrate for and continue
    // solving until we finish
    double elapsed_solve_time = 0;
    // this rate is updated at each time step and would typically vary with time
    state.SetCustomRateParameter("my photolysis rate", photo_rate);
    solver.CalculateRateConstants(state);

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
    }

    state.PrintState(time_step * (i + 1));
    photo_rate *= 1.5;
  }

  return 0;
}
