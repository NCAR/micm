#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/phase_transfer_process_builder.hpp>
#include <micm/process/process.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/process/rate_constant/user_defined_rate_constant.hpp>
#include <micm/process/transfer_coefficient/phase_transfer_coefficient.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <random>
using namespace micm;

template<class DenseMatrixPolicy>
void testProcessUpdateState(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 }, { "diffusion coefficient [m2 s-1]", 2.3e2 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  Phase gas_phase{ "gas", std::vector<micm::Species>{ foo, bar } };

  ArrheniusRateConstant rc1({ .A_ = 12.2, .C_ = 300.0 });
  SurfaceRateConstant rc2({ .label_ = "foo_surf", .species_ = foo });
  UserDefinedRateConstant rc3({ .label_ = "bar_user" });

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }
  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar'" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  DenseMatrixPolicy expected_rate_constants(number_of_grid_cells, 3, 0.0);
  std::vector<double> params = { 0.0, 0.0, 0.0 };
  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    state.conditions_[i_cell].temperature_ = get_double() * 285.0;
    state.conditions_[i_cell].pressure_ = get_double() * 101100.0;
    state.conditions_[i_cell].air_density_ = get_double() * 10.0;
    params[0] = get_double() * 1.0e-8;
    params[1] = get_double() * 1.0e5;
    params[2] = get_double() * 1.0e-2;
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["foo_surf.effective radius [m]"]] = params[0];
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["foo_surf.particle number concentration [# m-3]"]] = params[1];
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["bar_user"]] = params[2];
    std::vector<double>::const_iterator param_iter = params.begin();
    expected_rate_constants[i_cell][0] = rc1.Calculate(state.conditions_[i_cell], param_iter) * (state.conditions_[i_cell].air_density_ * 0.82);
    param_iter += rc1.SizeCustomParameters();
    expected_rate_constants[i_cell][1] = rc2.Calculate(state.conditions_[i_cell], param_iter);
    param_iter += rc2.SizeCustomParameters();
    expected_rate_constants[i_cell][2] =
        rc3.Calculate(state.conditions_[i_cell], param_iter) * (state.conditions_[i_cell].air_density_ * 0.82);
    param_iter += rc3.SizeCustomParameters();
  }

  Process::CalculateRateConstants(processes, state);

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
      EXPECT_EQ(state.rate_constants_[i_cell][i_rxn], expected_rate_constants[i_cell][i_rxn])
          << "grid cell " << i_cell << "; reaction " << i_rxn;
}

template<class T>
using Group1VectorMatrix = VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = VectorMatrix<T, 4>;

TEST(Process, Matrix)
{
  testProcessUpdateState<Matrix<double>>(5);
}

TEST(Process, VectorMatrix)
{
  testProcessUpdateState<Group1VectorMatrix<double>>(5);
  testProcessUpdateState<Group2VectorMatrix<double>>(5);
  testProcessUpdateState<Group3VectorMatrix<double>>(5);
  testProcessUpdateState<Group4VectorMatrix<double>>(5);
}

TEST(Process, BuildsChemicalReactionAndPhaseTransferProcess)
{
  auto O3 = Species("O3");
  auto NO = Species("NO");
  auto NO2 = Species("NO2");
  auto O2 = Species("O2");
  auto CO2 = Species{ "CO2" };
  auto H2O = Species{ "H2O" };
  auto Hplus = Species{ "H+" };
  auto CO32minus = Species{ "CO32-" };

  Phase gas_phase{ "gas", std::vector<Species>{ O3, NO, NO2, O2 } };
  Phase aqueous_phase{ "aqueous", std::vector<Species>{ CO2, H2O, Hplus, CO32minus } };

  // Build a ChemicalReaction
  Process chemical_reaction = ChemicalReactionBuilder()
                                  .SetReactants({ Species("O3"), Species("NO") })
                                  .SetProducts({ Yield(Species("NO2"), 1.0), Yield(Species("O2"), 1.0) })
                                  .SetRateConstant(ArrheniusRateConstant())
                                  .SetPhase(gas_phase)
                                  .Build();

  // Build a PhaseTransferProcess
  Process phase_transfer = PhaseTransferProcessBuilder()
                               .SetGasSpecies(gas_phase, { CO2 })
                               .SetCondensedSpecies(aqueous_phase, { Yield(Hplus, 2.0), Yield(CO32minus) })
                               .SetSolvent(aqueous_phase, H2O)
                               .SetTransferCoefficient(PhaseTransferCoefficient())
                               .Build();

  // Check that the first process is a ChemicalReaction
  std::visit(
      [](auto&& value)
      {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, ChemicalReaction>)
        {
          EXPECT_EQ(value.reactants_.size(), 2);
          EXPECT_EQ(value.products_[0].species_.name_, "NO2");
          EXPECT_EQ(value.phase_.name_, "gas");
        }
        else
        {
          FAIL() << "Expected ChemicalReaction, got different type";
        }
      },
      chemical_reaction.process_);

  // Check that the second process is a PhaseTransferProcess
  std::visit(
      [](auto&& value)
      {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, PhaseTransferProcess>)
        {
          EXPECT_EQ(value.gas_phase_.name_, "gas");
          EXPECT_EQ(value.condensed_phase_.name_, "aqueous");
          EXPECT_EQ(value.solvent_phase_.name_, "aqueous");
          EXPECT_EQ(value.gas_species_.size(), 1);
          EXPECT_EQ(value.condensed_species_.size(), 2);
          EXPECT_EQ(value.solvent_.name_, "H2O");
          EXPECT_EQ(value.condensed_species_[0].coefficient_, 2.0);
        }
        else
        {
          FAIL() << "Expected PhaseTransferProcess, got different type";
        }
      },
      phase_transfer.process_);
}

TEST(Process, ChemicalReactionCopyAssignmentSucceeds)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 }, { "diffusion coefficient [m2 s-1]", 2.3e2 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  Phase gas_phase{ "gas", std::vector<micm::Species>{ foo, bar } };
  ArrheniusRateConstant rc1({ .A_ = 12.2, .C_ = 300.0 });

  Process reaction = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1).SetPhase(gas_phase).Build();

  // Assign original to copy
  Process copy_reaction = reaction;

  std::visit(
      [](auto&& copy, auto&& original)
      {
        using T = std::decay_t<decltype(copy)>;
        using U = std::decay_t<decltype(original)>;
        if constexpr (std::is_same_v<T, ChemicalReaction> && std::is_same_v<U, ChemicalReaction>)
        {
          EXPECT_EQ(copy.reactants_[0].name_, original.reactants_[0].name_);
          EXPECT_EQ(copy.products_.size(), original.products_.size());
          EXPECT_EQ(copy.phase_.name_, original.phase_.name_);
          EXPECT_NE(copy.rate_constant_.get(), original.rate_constant_.get());
        }
        else
        {
          FAIL() << "Expected both variants to hold ChemicalReaction";
        }
      },
      copy_reaction.process_,
      reaction.process_);
}

TEST(Process, PhaseTransferProcessCopyAssignmentSucceeds)
{
  auto O3 = Species("O3");
  auto NO = Species("NO");
  auto NO2 = Species("NO2");
  auto O2 = Species("O2");
  auto CO2 = Species{ "CO2" };
  auto H2O = Species{ "H2O" };
  auto Hplus = Species{ "H+" };
  auto CO32minus = Species{ "CO32-" };

  Phase gas_phase{ "gas", std::vector<Species>{ O3, NO, NO2, O2 } };
  Phase aqueous_phase{ "aqueous", std::vector<Species>{ CO2, H2O, Hplus, CO32minus } };

  // Build a PhaseTransferProcess
  Process phase_transfer = PhaseTransferProcessBuilder()
                               .SetGasSpecies(gas_phase, { CO2 })
                               .SetCondensedSpecies(aqueous_phase, { Yield(Hplus, 2.0), Yield(CO32minus) })
                               .SetSolvent(aqueous_phase, H2O)
                               .SetTransferCoefficient(PhaseTransferCoefficient())
                               .Build();

  // Assign original to copy
  Process copy_process = phase_transfer;

  std::visit(
      [](auto&& copy, auto&& original)
      {
        using T = std::decay_t<decltype(copy)>;
        using U = std::decay_t<decltype(original)>;
        if constexpr (std::is_same_v<T, PhaseTransferProcess> && std::is_same_v<U, PhaseTransferProcess>)
        {
          EXPECT_EQ(copy.gas_phase_.name_, original.gas_phase_.name_);
          EXPECT_EQ(copy.condensed_phase_.name_, original.condensed_phase_.name_);
          EXPECT_EQ(copy.solvent_phase_.name_, original.solvent_phase_.name_);
          EXPECT_EQ(copy.gas_species_[0].name_, original.gas_species_[0].name_);
          EXPECT_EQ(copy.condensed_species_[1].species_.name_, original.condensed_species_[1].species_.name_);
          EXPECT_EQ(copy.solvent_.name_, original.solvent_.name_);
          EXPECT_NE(copy.coefficient_.get(), original.coefficient_.get());
        }
        else
        {
          FAIL() << "Expected both variants to hold PhaseTransferProcess";
        }
      },
      copy_process.process_,
      phase_transfer.process_);
}

TEST(Process, SurfaceRateConstantOnlyHasOneReactant)
{
  Species c("c", { { "molecular weight [kg mol-1]", 0.025 }, { "diffusion coefficient [m2 s-1]", 2.3e2 } });
  Species e("e");
  Phase gas_phase{ "gas", std::vector<micm::Species>{ c, e } };

  EXPECT_ANY_THROW(
      Process r = ChemicalReactionBuilder()
                      .SetReactants({ c, c })
                      .SetProducts({ Yield(e, 1) })
                      .SetRateConstant(SurfaceRateConstant({ .label_ = "c", .species_ = c, .reaction_probability_ = 0.90 }))
                      .SetPhase(gas_phase)
                      .Build(););
}
