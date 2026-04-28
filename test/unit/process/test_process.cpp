// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/phase_transfer_process_builder.hpp>
#include <micm/process/process.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/process/rate_constant/user_defined_rate_constant.hpp>
#include <micm/process/reaction_rate_store.hpp>
#include <micm/process/transfer_coefficient/henrys_law_constant.hpp>
#include <micm/util/constants.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#define _USE_MATH_DEFINES
#include <math.h>
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

using namespace micm;

/// @brief Sort processes in the same order that SolverBuilder uses, so rc_indices match.
static void SortLikeBuilder(std::vector<Process>& procs)
{
  std::stable_sort(
      procs.begin(),
      procs.end(),
      [](const Process& a, const Process& b)
      {
        auto order = [](const Process& p) -> int
        {
          const ChemicalReaction* rxn = std::get_if<ChemicalReaction>(&p.process_);
          if (!rxn)
            return 10;
          return std::visit(
              [](const auto& v) -> int
              {
                using T = std::decay_t<decltype(v)>;
                if constexpr (std::is_same_v<T, ArrheniusRateConstantParameters>)
                  return 0;
                else if constexpr (std::is_same_v<T, TroeRateConstantParameters>)
                  return 1;
                else if constexpr (std::is_same_v<T, TernaryChemicalActivationRateConstantParameters>)
                  return 2;
                else if constexpr (std::is_same_v<T, BranchedRateConstantParameters>)
                  return 3;
                else if constexpr (std::is_same_v<T, TunnelingRateConstantParameters>)
                  return 4;
                else if constexpr (std::is_same_v<T, TaylorSeriesRateConstantParameters>)
                  return 5;
                else if constexpr (std::is_same_v<T, ReversibleRateConstantParameters>)
                  return 6;
                else if constexpr (std::is_same_v<T, UserDefinedRateConstantParameters>)
                  return 7;
                else if constexpr (std::is_same_v<T, SurfaceRateConstantParameters>)
                  return 8;
                else
                  return 9;
              },
              rxn->rate_constant_);
        };
        return order(a) < order(b);
      });
}

/// @brief Return the custom parameter labels that BuildFrom assigns for a process.
static std::vector<std::string> CustomParamLabels(const Process& proc)
{
  if (const auto* rxn = std::get_if<ChemicalReaction>(&proc.process_))
  {
    if (const auto* p = std::get_if<UserDefinedRateConstantParameters>(&rxn->rate_constant_))
      return { p->label_ };
    if (const auto* p = std::get_if<SurfaceRateConstantParameters>(&rxn->rate_constant_))
      return { p->label_ + ".effective radius [m]",
               p->label_ + ".particle number concentration [# m-3]" };
  }
  return {};
}

template<class DenseMatrixPolicy>
void testProcessUpdateState(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  double foo_diff_coeff = 2.3e2;
  PhaseSpecies gas_foo(foo, foo_diff_coeff);
  PhaseSpecies gas_bar(bar);
  Phase gas_phase{ "gas", { gas_foo, gas_bar } };

  ArrheniusRateConstantParameters rc1_params{ .A_ = 12.2, .C_ = 300.0 };
  SurfaceRateConstantParameters   rc2_params{ .label_ = "foo_surf", .phase_species_ = gas_foo };
  UserDefinedRateConstantParameters rc3_params{ .label_ = "bar_user" };

  Process r1 = ChemicalReactionBuilder().SetReactants({ foo, bar }).SetRateConstant(rc1_params).SetPhase(gas_phase).Build();
  Process r2 = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(rc2_params).SetPhase(gas_phase).Build();
  Process r3 = ChemicalReactionBuilder().SetReactants({ bar }).SetRateConstant(rc3_params).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  // Sort as SolverBuilder does: Arrhenius(0), UserDefined(7), Surface(8)
  SortLikeBuilder(processes);

  // Build custom param labels in sorted order
  std::vector<std::string> param_labels;
  for (const auto& proc : processes)
  {
    auto labels = CustomParamLabels(proc);
    param_labels.insert(param_labels.end(), labels.begin(), labels.end());
  }

  // Build the ReactionRateConstantStore
  auto store = ReactionRateConstantStore::BuildFrom(processes);

  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    state.conditions_[i_cell].temperature_ = get_double() * 285.0;
    state.conditions_[i_cell].pressure_    = get_double() * 101100.0;
    state.conditions_[i_cell].air_density_ = get_double() * 10.0;

    double user_rate = get_double() * 1.0e-2;
    double radius    = get_double() * 1.0e-8;
    double num_conc  = get_double() * 1.0e5;

    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["bar_user"]]                               = user_rate;
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["foo_surf.effective radius [m]"]]          = radius;
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["foo_surf.particle number concentration [# m-3]"]] =
        num_conc;
  }

  ReactionRateConstantStore::EvaluateCpuRateConstants(store, state);
  ReactionRateConstantStore::CpuCalculateRateConstants(store, state);

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    const auto& cond     = state.conditions_[i_cell];
    double      user_rate = state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["bar_user"]];
    double      radius    = state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["foo_surf.effective radius [m]"]];
    double      num_conc  = state.custom_rate_parameters_[i_cell]
                                [state.custom_rate_parameter_map_["foo_surf.particle number concentration [# m-3]"]];

    // r1 (Arrhenius) at rc_index 0; bar is parameterized: air_density * 0.82
    double expected_arr = CalculateArrhenius(rc1_params, cond.temperature_, cond.pressure_);
    expected_arr *= (cond.air_density_ * 0.82);
    EXPECT_NEAR(state.rate_constants_[i_cell][0], expected_arr, 1.0e-12 * expected_arr)
        << "grid cell " << i_cell << "; Arrhenius reaction";

    // r3 (UserDefined) at rc_index = user_defined_offset; bar is parameterized
    double expected_ud = user_rate * (cond.air_density_ * 0.82);
    EXPECT_NEAR(state.rate_constants_[i_cell][store.user_defined_offset()], expected_ud, 1.0e-12 * expected_ud)
        << "grid cell " << i_cell << "; UserDefined reaction";

    // r2 (Surface) at rc_index = surface_offset; foo is not parameterized
    double mean_free_speed = std::sqrt(8.0 * constants::GAS_CONSTANT / (M_PI * 0.025) * cond.temperature_);
    double expected_surf   = 4.0 * num_conc * M_PI * radius * radius /
                             (radius / foo_diff_coeff + 4.0 / mean_free_speed);
    EXPECT_NEAR(state.rate_constants_[i_cell][store.surface_offset()], expected_surf, 1.0e-10 * expected_surf)
        << "grid cell " << i_cell << "; Surface reaction";
  }
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
  auto H2OCO3 = Species{ "H2CO3" };

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ O3, NO, NO2, O2 } };
  Phase aqueous_phase{ "aqueous", std::vector<PhaseSpecies>{ CO2, H2O, Hplus, CO32minus, H2OCO3 } };

  // Build a ChemicalReaction
  Process chemical_reaction = ChemicalReactionBuilder()
                                  .SetReactants({ Species("O3"), Species("NO") })
                                  .SetProducts({ StoichSpecies(Species("NO2"), 1.0), StoichSpecies(Species("O2"), 1.0) })
                                  .SetRateConstant(ArrheniusRateConstantParameters{})
                                  .SetPhase(gas_phase)
                                  .Build();

  // Build a PhaseTransferProcess
  Process phase_transfer = PhaseTransferProcessBuilder()
                               .SetGasSpecies(gas_phase, CO2)
                               .SetCondensedSpecies(aqueous_phase, H2OCO3)
                               .SetSolvent(aqueous_phase, H2O)
                               .SetTransferCoefficient(HenrysLawConstant())
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
          EXPECT_EQ(value.gas_species_.name_, "CO2");
          EXPECT_EQ(value.condensed_species_.name_, "H2CO3");
          EXPECT_EQ(value.solvent_.name_, "H2O");
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

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ foo, bar } };

  Process reaction = ChemicalReactionBuilder()
                         .SetReactants({ foo, bar })
                         .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 12.2, .C_ = 300.0 })
                         .SetPhase(gas_phase)
                         .Build();

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
          // With value semantics, rate_constant_ is copied by value — verify parameters match
          const auto* copy_arr = std::get_if<ArrheniusRateConstantParameters>(&copy.rate_constant_);
          const auto* orig_arr = std::get_if<ArrheniusRateConstantParameters>(&original.rate_constant_);
          EXPECT_NE(copy_arr, nullptr);
          EXPECT_NE(orig_arr, nullptr);
          if (copy_arr && orig_arr)
            EXPECT_EQ(copy_arr->A_, orig_arr->A_);
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
  auto H2OCO3 = Species{ "H2CO3" };

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ O3, NO, NO2, O2 } };
  Phase aqueous_phase{ "aqueous", std::vector<PhaseSpecies>{ CO2, H2O, Hplus, CO32minus, H2OCO3 } };

  // Build a PhaseTransferProcess
  Process phase_transfer = PhaseTransferProcessBuilder()
                               .SetGasSpecies(gas_phase, CO2)
                               .SetCondensedSpecies(aqueous_phase, H2OCO3)
                               .SetSolvent(aqueous_phase, H2O)
                               .SetTransferCoefficient(HenrysLawConstant())
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
          EXPECT_EQ(copy.gas_species_.name_, original.gas_species_.name_);
          EXPECT_EQ(copy.condensed_species_.name_, original.condensed_species_.name_);
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
