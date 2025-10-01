#include <micm/CPU.hpp>

#include <iomanip>
#include <iostream>
#include <map>

// Use our namespace so that this example is easier to read
using namespace micm;

// Conversion factor from moles m-3 to molecules cm-3 for consistency
// with the configuraion file
constexpr double MOLES_M3_TO_MOLECULES_CM3 = 1.0e-6 * constants::AVOGADRO_CONSTANT;

int main(const int argc, const char* argv[])
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species(
      "C",
      std::map<std::string, double>{ { "molecular weight [kg mol-1]", 0.025 }});
  auto d = Species("D");
  auto e = Species("E");
  auto f = Species("F");
  auto g = Species("G");
  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ a, b, c, d, e, f, g } };

  auto& phase_species_list = gas_phase.phase_species_;
  auto it = std::find_if(phase_species_list.begin(), phase_species_list.end(),
    [&c](const PhaseSpecies& ps) {
      return ps.species_.name_ == c.name_; });

  if (it == phase_species_list.end())
  {
    std::cout << "Species not found\n";
    return 1; // Failure
  }

  double c_diffusion_coefficient = 2.3e2;
  size_t surface_c_index = std::distance(phase_species_list.begin(), it);
  phase_species_list[surface_c_index].SetDiffusionCoefficient(c_diffusion_coefficient);

  Process r1 = ChemicalReactionBuilder()
                   .SetReactants({ a })
                   .SetProducts({ Yield(b, 1) })
                   .SetRateConstant(ArrheniusRateConstant({ .A_ = 2.15e-4, .B_ = 0, .C_ = 110 }))
                   .SetPhase(gas_phase)
                   .Build();

  // a branched reaction has two output pathways
  // this is represnted internal to micm as two different reactions
  auto branched_params = BranchedRateConstantParameters{ .X_ = 1.2, .Y_ = 204.3, .a0_ = 1.0e-3, .n_ = 2 };
  branched_params.branch_ = BranchedRateConstantParameters::Branch::Alkoxy;

  Process r2 = ChemicalReactionBuilder()
                   .SetReactants({ b })
                   .SetProducts({ Yield(c, 1) })
                   .SetRateConstant(BranchedRateConstant(branched_params))
                   .SetPhase(gas_phase)
                   .Build();

  branched_params.branch_ = BranchedRateConstantParameters::Branch::Nitrate;
  Process r3 = ChemicalReactionBuilder()
                   .SetReactants({ b })
                   .SetProducts({ Yield(d, 1) })
                   .SetRateConstant(BranchedRateConstant(branched_params))
                   .SetPhase(gas_phase)
                   .Build();

  // A surface rate constant also needs to know the effective radius and particle number concentration
  // we will set those later
  Process r4 = ChemicalReactionBuilder()
                   .SetReactants({ c })
                   .SetProducts({ Yield(e, 1) })
                   .SetRateConstant(SurfaceRateConstant({ .label_ = "C", .phase_species_ = phase_species_list[surface_c_index], .reaction_probability_ = 0.90 }))
                   .SetPhase(gas_phase)
                   .Build();

  Process r5 = ChemicalReactionBuilder()
                   .SetReactants({ d })
                   .SetProducts({ Yield(f, 2) })
                   .SetRateConstant(TernaryChemicalActivationRateConstant({ .k0_A_ = 1.2,
                                                                            .k0_B_ = 2.3,
                                                                            .k0_C_ = 302.3,
                                                                            .kinf_A_ = 2.6 / MOLES_M3_TO_MOLECULES_CM3,
                                                                            .kinf_B_ = -3.1,
                                                                            .kinf_C_ = 402.1,
                                                                            .Fc_ = 0.9,
                                                                            .N_ = 1.2 }))
                   .SetPhase(gas_phase)
                   .Build();

  // to have a stoichiemetric coefficient of more than one for reactants,
  // list the reactant that many times
  Process r6 =
      ChemicalReactionBuilder()
          .SetReactants({ e, e })
          .SetProducts({ Yield(g, 1) })
          .SetRateConstant(TroeRateConstant({ .k0_A_ = 1.2e4 * MOLES_M3_TO_MOLECULES_CM3 * MOLES_M3_TO_MOLECULES_CM3,
                                              .k0_B_ = 167.0,
                                              .k0_C_ = 3.0,
                                              .kinf_A_ = 136.0 * MOLES_M3_TO_MOLECULES_CM3,
                                              .kinf_B_ = 5.0,
                                              .kinf_C_ = 24.0,
                                              .Fc_ = 0.9,
                                              .N_ = 0.8 }))
          .SetPhase(gas_phase)
          .Build();

  Process r7 = ChemicalReactionBuilder()
                   .SetReactants({ f })
                   .SetProducts({ Yield(g, 1) })
                   .SetRateConstant(TunnelingRateConstant({ .A_ = 1.2, .B_ = 2.3, .C_ = 302.3 }))
                   .SetPhase(gas_phase)
                   .Build();

  Process r8 = ChemicalReactionBuilder()
                   .SetReactants({ c })
                   .SetProducts({ Yield(g, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "my photolysis rate" }))
                   .SetPhase(gas_phase)
                   .Build();

  Process r9 = ChemicalReactionBuilder()
                   .SetProducts({ Yield(a, 1) })
                   .SetRateConstant(UserDefinedRateConstant({ .label_ = "my emission rate" }))
                   .SetPhase(gas_phase)
                   .Build();

  Process r10 = ChemicalReactionBuilder()
                    .SetReactants({ b })
                    .SetRateConstant(UserDefinedRateConstant({ .label_ = "my loss rate" }))
                    .SetPhase(gas_phase)
                    .Build();

  auto chemical_system = System(micm::SystemParameters{ .gas_phase_ = gas_phase });
  auto reactions = std::vector<micm::Process>{ r1, r2, r3, r4, r5, r6, r7, r8, r9, r10 };

  auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(
                    micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(chemical_system)
                    .SetReactions(reactions)
                    .Build();
  State state = solver.GetState();

  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa
  state.conditions_[0].CalculateIdealAirDensity();

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
      elapsed_solve_time += result.final_time_;
    }

    state.PrintState(time_step * (i + 1));
    photo_rate *= 1.5;
  }

  return 0;
}
