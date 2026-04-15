// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/branched_rate_constant.hpp>
#include <micm/process/rate_constant/lambda_rate_constant.hpp>
#include <micm/process/rate_constant/reversible_rate_constant.hpp>
#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/process/rate_constant/taylor_series_rate_constant.hpp>
#include <micm/process/rate_constant/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/rate_constant/troe_rate_constant.hpp>
#include <micm/process/rate_constant/tunneling_rate_constant.hpp>
#include <micm/process/rate_constant/user_defined_rate_constant.hpp>
#include <micm/process/reaction_rate_store.hpp>
#include <micm/util/constants.hpp>

#include <gtest/gtest.h>

#define _USE_MATH_DEFINES
#include <math.h>
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

using namespace micm;

namespace
{
  // Helper to create a minimal gas phase with species
  Phase MakeGasPhase(std::vector<Species> species)
  {
    std::vector<PhaseSpecies> ps;
    ps.reserve(species.size());
    for (auto& s : species)
      ps.emplace_back(s);
    return Phase{ "gas", ps };
  }

  // Helper to create a gas phase with one species carrying a diffusion coefficient
  Phase MakeGasPhaseWithDiffusion(const Species& sp, double diffusion_coefficient)
  {
    PhaseSpecies ps(sp, diffusion_coefficient);
    return Phase{ "gas", { ps } };
  }

  // Helper: build processes list sorted exactly as SolverBuilder would
  // (reactions already sorted by RateConstantTypeOrder before calling BuildFrom)
  void SortByTypeOrder(std::vector<Process>& procs)
  {
    auto type_order = [](const Process& p) -> int
    {
      const auto* rxn = std::get_if<ChemicalReaction>(&p.process_);
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
    std::stable_sort(procs.begin(), procs.end(), [&](const Process& a, const Process& b)
                     { return type_order(a) < type_order(b); });
  }
}  // namespace

// ============================================================
// Offset arithmetic
// ============================================================

TEST(ReactionRateStore, OffsetsAreContiguousCumulativeSizes)
{
  Species a("a"), b("b"), c("c", { { "molecular weight [kg mol-1]", 0.025 } });
  double c_diff = 1.0e-5;
  PhaseSpecies gas_c(c, c_diff);
  Phase gas{ "gas", { PhaseSpecies(a), PhaseSpecies(b), gas_c } };

  std::vector<Process> procs;
  procs.push_back(ChemicalReactionBuilder()
                      .SetReactants({ a })
                      .SetProducts({ StoichSpecies(b, 1) })
                      .SetRateConstant(ArrheniusRateConstantParameters{})
                      .SetPhase(gas)
                      .Build());
  procs.push_back(ChemicalReactionBuilder()
                      .SetReactants({ a })
                      .SetProducts({ StoichSpecies(b, 1) })
                      .SetRateConstant(ArrheniusRateConstantParameters{})
                      .SetPhase(gas)
                      .Build());
  procs.push_back(ChemicalReactionBuilder()
                      .SetReactants({ a })
                      .SetProducts({ StoichSpecies(b, 1) })
                      .SetRateConstant(TroeRateConstantParameters{ .k0_A_ = 1.0, .kinf_A_ = 1.0 })
                      .SetPhase(gas)
                      .Build());
  procs.push_back(ChemicalReactionBuilder()
                      .SetReactants({ a })
                      .SetProducts({ StoichSpecies(b, 1) })
                      .SetRateConstant(TunnelingRateConstantParameters{})
                      .SetPhase(gas)
                      .Build());
  procs.push_back(ChemicalReactionBuilder()
                      .SetReactants({ a })
                      .SetProducts({ StoichSpecies(b, 1) })
                      .SetRateConstant(UserDefinedRateConstantParameters{ .label_ = "ud1" })
                      .SetPhase(gas)
                      .Build());
  procs.push_back(ChemicalReactionBuilder()
                      .SetReactants({ c })
                      .SetProducts({ StoichSpecies(b, 1) })
                      .SetRateConstant(SurfaceRateConstantParameters{ .label_ = "surf1", .phase_species_ = gas_c })
                      .SetPhase(gas)
                      .Build());

  SortByTypeOrder(procs);
  auto store = ReactionRateStore::BuildFrom(procs);

  // Arrhenius: 2, Troe: 1, Ternary: 0, Branched: 0, Tunneling: 1, Taylor: 0, Reversible: 0, UserDefined: 1, Surface: 1
  EXPECT_EQ(store.arrhenius.size(), 2u);
  EXPECT_EQ(store.troe.size(), 1u);
  EXPECT_EQ(store.ternary.size(), 0u);
  EXPECT_EQ(store.branched.size(), 0u);
  EXPECT_EQ(store.tunneling.size(), 1u);
  EXPECT_EQ(store.taylor.size(), 0u);
  EXPECT_EQ(store.reversible.size(), 0u);
  EXPECT_EQ(store.user_defined.size(), 1u);
  EXPECT_EQ(store.surface.size(), 1u);

  // Verify cumulative offsets
  EXPECT_EQ(store.troe_offset(), 2u);
  EXPECT_EQ(store.ternary_offset(), 3u);
  EXPECT_EQ(store.branched_offset(), 3u);
  EXPECT_EQ(store.tunneling_offset(), 3u);
  EXPECT_EQ(store.taylor_offset(), 4u);
  EXPECT_EQ(store.reversible_offset(), 4u);
  EXPECT_EQ(store.user_defined_offset(), 4u);
  EXPECT_EQ(store.surface_offset(), 5u);
  EXPECT_EQ(store.lambda_offset(), 6u);

  // Total rate constants (lambda_offset + lambda_entries.size()) equals process count
  EXPECT_EQ(store.lambda_offset() + store.lambda_entries.size(), procs.size());
}

// ============================================================
// Single Arrhenius: parameters preserved
// ============================================================

TEST(ReactionRateStore, ArrheniusParametersPreserved)
{
  Species a("a"), b("b");
  Phase gas = MakeGasPhase({ a, b });

  ArrheniusRateConstantParameters params{ .A_ = 2.15e-4, .B_ = 1.2, .C_ = 110.0, .D_ = 300.0, .E_ = 0.0 };
  std::vector<Process> procs{ ChemicalReactionBuilder()
                                   .SetReactants({ a })
                                   .SetProducts({ StoichSpecies(b, 1) })
                                   .SetRateConstant(params)
                                   .SetPhase(gas)
                                   .Build() };

  auto store = ReactionRateStore::BuildFrom(procs);
  ASSERT_EQ(store.arrhenius.size(), 1u);
  EXPECT_DOUBLE_EQ(store.arrhenius[0].A_, 2.15e-4);
  EXPECT_DOUBLE_EQ(store.arrhenius[0].B_, 1.2);
  EXPECT_DOUBLE_EQ(store.arrhenius[0].C_, 110.0);
  EXPECT_DOUBLE_EQ(store.arrhenius[0].D_, 300.0);
}

// ============================================================
// Branched: derived fields k0_ and z_ are computed
// ============================================================

TEST(ReactionRateStore, BranchedDerivedFieldsComputed)
{
  Species a("a"), b("b");
  Phase gas = MakeGasPhase({ a, b });

  BranchedRateConstantParameters params{
    .branch_ = BranchedRateConstantParameters::Branch::Alkoxy, .X_ = 1.0, .Y_ = 0.0, .a0_ = 0.5, .n_ = 3
  };
  std::vector<Process> procs{ ChemicalReactionBuilder()
                                   .SetReactants({ a })
                                   .SetProducts({ StoichSpecies(b, 1) })
                                   .SetRateConstant(params)
                                   .SetPhase(gas)
                                   .Build() };

  auto store = ReactionRateStore::BuildFrom(procs);
  ASSERT_EQ(store.branched.size(), 1u);

  // k0_ = 2e-22 * N_A * 1e-6 * exp(n)
  double expected_k0 = 2.0e-22 * constants::AVOGADRO_CONSTANT * 1.0e-6 * std::exp(3.0);
  EXPECT_NEAR(store.branched[0].k0_, expected_k0, 1.0e-10 * expected_k0);

  // z_ = A_val * (1 - a0) / a0
  double air_ref = 2.45e19 / constants::AVOGADRO_CONSTANT * 1.0e6;
  double a_val = expected_k0 * air_ref;
  double b_val = 0.43 * std::pow(293.0 / 298.0, -8.0);
  double A_val = a_val / (1.0 + a_val / b_val) * std::pow(0.41, 1.0 / (1.0 + std::pow(std::log10(a_val / b_val), 2.0)));
  double expected_z = A_val * (1.0 - 0.5) / 0.5;
  EXPECT_NEAR(store.branched[0].z_, expected_z, 1.0e-10 * std::abs(expected_z));
}

// ============================================================
// UserDefined: custom_param_index_ assigned correctly
// ============================================================

TEST(ReactionRateStore, UserDefinedCustomParamIndex)
{
  Species a("a"), b("b");
  Phase gas = MakeGasPhase({ a, b });

  UserDefinedRateConstantParameters p1{ .label_ = "rate1", .scaling_factor_ = 2.0 };
  UserDefinedRateConstantParameters p2{ .label_ = "rate2", .scaling_factor_ = 0.5 };
  std::vector<Process> procs{
    ChemicalReactionBuilder().SetReactants({ a }).SetProducts({ StoichSpecies(b, 1) }).SetRateConstant(p1).SetPhase(gas).Build(),
    ChemicalReactionBuilder().SetReactants({ a }).SetProducts({ StoichSpecies(b, 1) }).SetRateConstant(p2).SetPhase(gas).Build()
  };

  auto store = ReactionRateStore::BuildFrom(procs);
  ASSERT_EQ(store.user_defined.size(), 2u);

  EXPECT_EQ(store.user_defined[0].custom_param_index_, 0u);
  EXPECT_DOUBLE_EQ(store.user_defined[0].scaling_factor_, 2.0);

  EXPECT_EQ(store.user_defined[1].custom_param_index_, 1u);
  EXPECT_DOUBLE_EQ(store.user_defined[1].scaling_factor_, 0.5);
}

// ============================================================
// Surface: data fields and custom_param_base_index_ assigned
// ============================================================

TEST(ReactionRateStore, SurfaceDataFieldsAndCustomParamIndex)
{
  double mw = 0.025;
  double diff_coeff = 2.3e2;
  double prob = 0.74;

  Species c("c", { { "molecular weight [kg mol-1]", mw } });
  Species b("b");
  PhaseSpecies gas_c(c, diff_coeff);
  Phase gas{ "gas", { gas_c, PhaseSpecies(b) } };

  SurfaceRateConstantParameters params{ .label_ = "surf", .phase_species_ = gas_c, .reaction_probability_ = prob };
  std::vector<Process> procs{ ChemicalReactionBuilder()
                                   .SetReactants({ c })
                                   .SetProducts({ StoichSpecies(b, 1) })
                                   .SetRateConstant(params)
                                   .SetPhase(gas)
                                   .Build() };

  auto store = ReactionRateStore::BuildFrom(procs);
  ASSERT_EQ(store.surface.size(), 1u);

  EXPECT_DOUBLE_EQ(store.surface[0].diffusion_coefficient_, diff_coeff);
  EXPECT_NEAR(store.surface[0].mean_free_speed_factor_, 8.0 * constants::GAS_CONSTANT / (M_PI * mw), 1.0e-14);
  EXPECT_DOUBLE_EQ(store.surface[0].reaction_probability_, prob);
  EXPECT_EQ(store.surface[0].custom_param_base_index_, 0u);
}

TEST(ReactionRateStore, SurfaceCustomParamIndexAfterUserDefined)
{
  double mw = 0.025;
  double diff_coeff = 1.0e-5;

  Species a("a"), b("b"), c("c", { { "molecular weight [kg mol-1]", mw } });
  PhaseSpecies gas_c(c, diff_coeff);
  Phase gas{ "gas", { PhaseSpecies(a), PhaseSpecies(b), gas_c } };

  UserDefinedRateConstantParameters ud_p{ .label_ = "ud" };
  SurfaceRateConstantParameters surf_p{ .label_ = "surf", .phase_species_ = gas_c };

  // Sorted: UserDefined (7) < Surface (8)
  std::vector<Process> procs{
    ChemicalReactionBuilder().SetReactants({ a }).SetProducts({ StoichSpecies(b, 1) }).SetRateConstant(ud_p).SetPhase(gas).Build(),
    ChemicalReactionBuilder().SetReactants({ c }).SetProducts({ StoichSpecies(b, 1) }).SetRateConstant(surf_p).SetPhase(gas).Build()
  };

  auto store = ReactionRateStore::BuildFrom(procs);
  ASSERT_EQ(store.user_defined.size(), 1u);
  ASSERT_EQ(store.surface.size(), 1u);

  // UserDefined takes slot 0; Surface takes slots 1 (radius) and 2 (num_conc)
  EXPECT_EQ(store.user_defined[0].custom_param_index_, 0u);
  EXPECT_EQ(store.surface[0].custom_param_base_index_, 1u);
}

// ============================================================
// Lambda: lambda_entries populated with correct rc_index
// ============================================================

TEST(ReactionRateStore, LambdaEntriesRcIndex)
{
  Species a("a"), b("b");
  Phase gas = MakeGasPhase({ a, b });

  ArrheniusRateConstantParameters arr{};
  LambdaRateConstantParameters lam{
    .label_ = "lam",
    .lambda_function_ = [](const Conditions& c) { return c.temperature_ * 1.0e-3; }
  };

  // Sorted: Arrhenius (0) before Lambda (9)
  std::vector<Process> procs{
    ChemicalReactionBuilder().SetReactants({ a }).SetProducts({ StoichSpecies(b, 1) }).SetRateConstant(arr).SetPhase(gas).Build(),
    ChemicalReactionBuilder().SetReactants({ a }).SetProducts({ StoichSpecies(b, 1) }).SetRateConstant(lam).SetPhase(gas).Build()
  };

  auto store = ReactionRateStore::BuildFrom(procs);

  EXPECT_EQ(store.arrhenius.size(), 1u);
  ASSERT_EQ(store.lambda_entries.size(), 1u);

  // Lambda is the second process → rc_index = 1
  EXPECT_EQ(store.lambda_entries[0].rc_index, 1u);
  // Pointer should be non-null and function should work
  ASSERT_NE(store.lambda_entries[0].source, nullptr);
  Conditions cond{ .temperature_ = 300.0 };
  EXPECT_NEAR(store.lambda_entries[0].source->lambda_function_(cond), 0.3, 1.0e-14);
}

// ============================================================
// Parameterized multipliers: one per parameterized reactant reaction
// ============================================================

TEST(ReactionRateStore, ParameterizedMultipliers)
{
  Species a("a"), b("b");
  b.parameterize_ = [](const Conditions& c) { return c.air_density_ * 2.0; };
  Phase gas = MakeGasPhase({ a, b });

  ArrheniusRateConstantParameters arr_a{};  // no parameterized reactants
  ArrheniusRateConstantParameters arr_b{};  // b is parameterized

  std::vector<Process> procs{
    ChemicalReactionBuilder().SetReactants({ a }).SetProducts({ StoichSpecies(b, 1) }).SetRateConstant(arr_a).SetPhase(gas).Build(),
    ChemicalReactionBuilder().SetReactants({ b }).SetProducts({ StoichSpecies(a, 1) }).SetRateConstant(arr_b).SetPhase(gas).Build()
  };

  auto store = ReactionRateStore::BuildFrom(procs);

  // Only the second reaction has a parameterized reactant
  ASSERT_EQ(store.parameterized_multipliers.size(), 1u);
  EXPECT_EQ(store.parameterized_multipliers[0].rc_index, 1u);

  Conditions cond{ .air_density_ = 5.0 };
  EXPECT_NEAR(store.parameterized_multipliers[0].evaluate(cond), 10.0, 1.0e-14);
}

// ============================================================
// Total rate constants matches process count for mixed mechanism
// ============================================================

TEST(ReactionRateStore, TotalRateConstantsMatchesProcessCount)
{
  Species a("a"), b("b"), c("c", { { "molecular weight [kg mol-1]", 0.025 } });
  PhaseSpecies gas_c(c, 1.0e-5);
  Phase gas{ "gas", { PhaseSpecies(a), PhaseSpecies(b), gas_c } };

  LambdaRateConstantParameters lam{
    .label_ = "lam", .lambda_function_ = [](const Conditions&) { return 1.0; }
  };
  SurfaceRateConstantParameters surf{ .label_ = "surf", .phase_species_ = gas_c };

  std::vector<Process> procs;
  procs.push_back(ChemicalReactionBuilder()
                      .SetReactants({ a })
                      .SetProducts({ StoichSpecies(b, 1) })
                      .SetRateConstant(ArrheniusRateConstantParameters{})
                      .SetPhase(gas)
                      .Build());
  procs.push_back(ChemicalReactionBuilder()
                      .SetReactants({ a })
                      .SetProducts({ StoichSpecies(b, 1) })
                      .SetRateConstant(TroeRateConstantParameters{ .k0_A_ = 1.0, .kinf_A_ = 1.0 })
                      .SetPhase(gas)
                      .Build());
  procs.push_back(ChemicalReactionBuilder()
                      .SetReactants({ a })
                      .SetProducts({ StoichSpecies(b, 1) })
                      .SetRateConstant(UserDefinedRateConstantParameters{ .label_ = "ud" })
                      .SetPhase(gas)
                      .Build());
  procs.push_back(ChemicalReactionBuilder()
                      .SetReactants({ c })
                      .SetProducts({ StoichSpecies(b, 1) })
                      .SetRateConstant(surf)
                      .SetPhase(gas)
                      .Build());
  procs.push_back(ChemicalReactionBuilder()
                      .SetReactants({ a })
                      .SetProducts({ StoichSpecies(b, 1) })
                      .SetRateConstant(lam)
                      .SetPhase(gas)
                      .Build());

  SortByTypeOrder(procs);
  auto store = ReactionRateStore::BuildFrom(procs);

  // Total = lambda_offset() + lambda_entries.size() = number of chemical reactions
  std::size_t n_rxn = 0;
  for (const auto& p : procs)
    if (std::holds_alternative<ChemicalReaction>(p.process_))
      ++n_rxn;

  EXPECT_EQ(store.lambda_offset() + store.lambda_entries.size(), n_rxn);
}

// ============================================================
// Error: missing diffusion coefficient
// ============================================================

TEST(ReactionRateStore, SurfaceMissingDiffusionCoefficientThrows)
{
  Species c("c", { { "molecular weight [kg mol-1]", 0.025 } });
  Species b("b");
  PhaseSpecies gas_c_no_diff(c);  // no diffusion coefficient
  Phase gas{ "gas", { gas_c_no_diff, PhaseSpecies(b) } };

  SurfaceRateConstantParameters params{ .label_ = "surf", .phase_species_ = gas_c_no_diff };
  std::vector<Process> procs{ ChemicalReactionBuilder()
                                   .SetReactants({ c })
                                   .SetProducts({ StoichSpecies(b, 1) })
                                   .SetRateConstant(params)
                                   .SetPhase(gas)
                                   .Build() };

  EXPECT_THROW({ auto store = ReactionRateStore::BuildFrom(procs); }, MicmException);
}

// ============================================================
// Error: missing molecular weight
// ============================================================

TEST(ReactionRateStore, SurfaceMissingMolecularWeightThrows)
{
  Species c("c");  // no properties at all
  Species b("b");
  PhaseSpecies gas_c(c, 1.0e-5);
  Phase gas{ "gas", { gas_c, PhaseSpecies(b) } };

  SurfaceRateConstantParameters params{ .label_ = "surf", .phase_species_ = gas_c };
  std::vector<Process> procs{ ChemicalReactionBuilder()
                                   .SetReactants({ c })
                                   .SetProducts({ StoichSpecies(b, 1) })
                                   .SetRateConstant(params)
                                   .SetPhase(gas)
                                   .Build() };

  EXPECT_THROW({ auto store = ReactionRateStore::BuildFrom(procs); }, MicmException);
}
