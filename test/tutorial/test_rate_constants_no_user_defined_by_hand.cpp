#include <iomanip>
#include <iostream>
#include <map>

// Each rate constant is in its own header file
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

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

  Process r1 = Process::create()
                   .reactants({ a })
                   .products({ yields(b, 1) })
                   .rate_constant(ArrheniusRateConstant({ .A_ = 2.15e-1, .B_ = 0, .C_ = 110 }))
                   .phase(gas_phase);

  // a branched reaction has two output pathways
  // this is represnted internal to micm as two different reactions
  auto branched_params = BranchedRateConstantParameters{ .X_ = 1.2, .Y_ = 204.3, .a0_ = 1.0e-3, .n_ = 2 };
  branched_params.branch_ = BranchedRateConstantParameters::Branch::Alkoxy;

  Process r2 = Process::create()
                   .reactants({ b })
                   .products({ yields(c, 1) })
                   .rate_constant(BranchedRateConstant(branched_params))
                   .phase(gas_phase);

  branched_params.branch_ = BranchedRateConstantParameters::Branch::Nitrate;
  Process r3 = Process::create()
                   .reactants({ b })
                   .products({ yields(d, 1) })
                   .rate_constant(BranchedRateConstant(branched_params))
                   .phase(gas_phase);

  // A surface rate constant also needs to know the effective radius and particle number concentration
  // we will set those later
  Process r4 = Process::create()
                   .reactants({ c })
                   .products({ yields(e, 1) })
                   .rate_constant(SurfaceRateConstant({ .label_ = "C", .species_ = c, .reaction_probability_ = 0.90 }))
                   .phase(gas_phase);

  Process r5 = Process::create()
                   .reactants({ d })
                   .products({ yields(f, 2) })
                   .rate_constant(TernaryChemicalActivationRateConstant({ .k0_A_ = 1.2,
                                                                          .k0_B_ = 2.3,
                                                                          .k0_C_ = 302.3,
                                                                          .kinf_A_ = 2.6,
                                                                          .kinf_B_ = -3.1,
                                                                          .kinf_C_ = 402.1,
                                                                          .Fc_ = 0.9,
                                                                          .N_ = 1.2 }))
                   .phase(gas_phase);

  // to have a stoichiemetric coefficient of more than one for reactants,
  // list the reactant that many times
  Process r6 = Process::create()
                   .reactants({ e, e })
                   .products({ yields(g, 1) })
                   .rate_constant(TroeRateConstant({ .k0_A_ = 1.2e4,
                                                     .k0_B_ = 167.0,
                                                     .k0_C_ = 3.0,
                                                     .kinf_A_ = 136.0,
                                                     .kinf_B_ = 5.0,
                                                     .kinf_C_ = 24.0,
                                                     .Fc_ = 0.9,
                                                     .N_ = 0.8 }))
                   .phase(gas_phase);

  Process r7 = Process::create()
                   .reactants({ f })
                   .products({ yields(g, 1) })
                   .rate_constant(TunnelingRateConstant({ .A_ = 1.2, .B_ = 2.3, .C_ = 302.3 }))
                   .phase(gas_phase);

  auto chemical_system = System(micm::SystemParameters{ .gas_phase_ = gas_phase });
  auto reactions = std::vector<micm::Process>{ r1, r2, r3, r4, r5, r6, r7 };

  RosenbrockSolver<> solver{ chemical_system, reactions, RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };
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

  // choose a timestep and print the initial state
  double time_step = 500;  // s

  state.PrintHeader();
  state.PrintState(0);

  // solve for ten iterations
  for (int i = 0; i < 10; ++i)
  {
    // Depending on how stiff the system is
    // the solver integration step may not be able to solve for the full time step
    // so we need to track how much time the solver was able to integrate for and continue
    // solving until we finish
    double elapsed_solve_time = 0;

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
      // std::cout << "solver state: " << StateToString(result.state_) << std::endl;
      state.variables_ = result.result_;
    }

    state.PrintState(time_step * (i + 1));
  }

  return 0;
}
