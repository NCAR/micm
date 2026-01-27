// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/process.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>

#include <gtest/gtest.h>

using namespace micm;

TEST(ChemicalReactionBuilder, SetAerosolScopeAppliesCorrectScopingToReactants)
{
  // Create species
  auto CO2 = Species{ "CO2" };
  auto H2O = Species{ "H2O" };

  // Create aqueous phase
  Phase aqueous_phase{ "aqueous", std::vector<PhaseSpecies>{ CO2, H2O } };

  // Create a reaction with aerosol scoping
  auto rate_constant = ArrheniusRateConstant{ { .A_ = 1.0 } };

  Process reaction = ChemicalReactionBuilder()
                         .SetAerosolScope("accumulation", aqueous_phase)
                         .SetReactants({ CO2 })
                         .SetProducts({ StoichSpecies(H2O, 1.0) })
                         .SetRateConstant(rate_constant)
                         .Build();

  // Verify the reaction was created and reactant name was scoped
  auto* chem_reaction = std::get_if<ChemicalReaction>(&reaction.process_);
  ASSERT_NE(chem_reaction, nullptr);
  ASSERT_EQ(chem_reaction->reactants_.size(), 1);
  EXPECT_EQ(chem_reaction->reactants_[0].name_, "accumulation.aqueous.CO2");

  // Verify product name was scoped
  ASSERT_EQ(chem_reaction->products_.size(), 1);
  EXPECT_EQ(chem_reaction->products_[0].species_.name_, "accumulation.aqueous.H2O");
}

TEST(ChemicalReactionBuilder, SetAerosolScopeAppliesCorrectScopingToProducts)
{
  // Create species
  auto OH = Species{ "OH-" };
  auto Hplus = Species{ "H+" };
  auto H2O = Species{ "H2O" };

  // Create aqueous phase
  Phase aqueous_phase{ "aqueous", std::vector<PhaseSpecies>{ OH, Hplus, H2O } };

  // Create a reaction with aerosol scoping
  auto rate_constant = ArrheniusRateConstant{ { .A_ = 1.14e-2 } };

  Process reaction = ChemicalReactionBuilder()
                         .SetAerosolScope("aitken", aqueous_phase)
                         .SetReactants({ H2O })
                         .SetProducts({ StoichSpecies(OH, 1.0), StoichSpecies(Hplus, 1.0) })
                         .SetRateConstant(rate_constant)
                         .Build();

  // Verify the reaction was created
  auto* chem_reaction = std::get_if<ChemicalReaction>(&reaction.process_);
  ASSERT_NE(chem_reaction, nullptr);

  // Verify reactant name was scoped
  ASSERT_EQ(chem_reaction->reactants_.size(), 1);
  EXPECT_EQ(chem_reaction->reactants_[0].name_, "aitken.aqueous.H2O");

  // Verify product names were scoped
  ASSERT_EQ(chem_reaction->products_.size(), 2);
  EXPECT_EQ(chem_reaction->products_[0].species_.name_, "aitken.aqueous.OH-");
  EXPECT_EQ(chem_reaction->products_[1].species_.name_, "aitken.aqueous.H+");
}

TEST(ChemicalReactionBuilder, SetAerosolScopeWithMultipleReactantsAndProducts)
{
  // Create species
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };
  auto D = Species{ "D" };

  // Create phase
  Phase organic_phase{ "organic", std::vector<PhaseSpecies>{ A, B, C, D } };

  // Create a reaction with aerosol scoping
  auto rate_constant = ArrheniusRateConstant{ { .A_ = 2.5 } };

  Process reaction = ChemicalReactionBuilder()
                         .SetAerosolScope("coarse", organic_phase)
                         .SetReactants({ A, B })
                         .SetProducts({ StoichSpecies(C, 1.0), StoichSpecies(D, 2.0) })
                         .SetRateConstant(rate_constant)
                         .Build();

  // Verify the reaction was created
  auto* chem_reaction = std::get_if<ChemicalReaction>(&reaction.process_);
  ASSERT_NE(chem_reaction, nullptr);

  // Verify reactant names were scoped
  ASSERT_EQ(chem_reaction->reactants_.size(), 2);
  EXPECT_EQ(chem_reaction->reactants_[0].name_, "coarse.organic.A");
  EXPECT_EQ(chem_reaction->reactants_[1].name_, "coarse.organic.B");

  // Verify product names and yields were preserved
  ASSERT_EQ(chem_reaction->products_.size(), 2);
  EXPECT_EQ(chem_reaction->products_[0].species_.name_, "coarse.organic.C");
  EXPECT_EQ(chem_reaction->products_[0].coefficient_, 1.0);
  EXPECT_EQ(chem_reaction->products_[1].species_.name_, "coarse.organic.D");
  EXPECT_EQ(chem_reaction->products_[1].coefficient_, 2.0);
}

TEST(ChemicalReactionBuilder, SetAerosolScopeDifferentPhaseNames)
{
  // Test that different phase names create different scopes
  auto species = Species{ "ABC" };

  Phase phase1{ "aqueous", std::vector<PhaseSpecies>{ species } };
  Phase phase2{ "organic", std::vector<PhaseSpecies>{ species } };

  auto rate_constant = ArrheniusRateConstant{ { .A_ = 1.0 } };

  // Reaction with phase1
  Process reaction1 = ChemicalReactionBuilder()
                          .SetAerosolScope("aitken", phase1)
                          .SetReactants({ species })
                          .SetProducts({})
                          .SetRateConstant(rate_constant)
                          .Build();

  // Reaction with phase2
  Process reaction2 = ChemicalReactionBuilder()
                          .SetAerosolScope("dust", phase2)
                          .SetReactants({ species })
                          .SetProducts({})
                          .SetRateConstant(rate_constant)
                          .Build();

  auto* chem_reaction1 = std::get_if<ChemicalReaction>(&reaction1.process_);
  auto* chem_reaction2 = std::get_if<ChemicalReaction>(&reaction2.process_);

  EXPECT_EQ(chem_reaction1->reactants_[0].name_, "aitken.aqueous.ABC");
  EXPECT_EQ(chem_reaction2->reactants_[0].name_, "dust.organic.ABC");
}

// ============================================================================
// Tests for mutual exclusivity of SetPhase and SetAerosolScope
// ============================================================================

TEST(ChemicalReactionBuilder, SetPhaseAfterSetAerosolScopeThrowsError)
{
  auto CO2 = Species{ "CO2" };
  Phase aqueous_phase{ "aqueous", std::vector<PhaseSpecies>{ CO2 } };
  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ CO2 } };
  auto rate_constant = ArrheniusRateConstant{ { .A_ = 1.0 } };

  EXPECT_THROW(
      {
        ChemicalReactionBuilder()
            .SetAerosolScope("accumulation", aqueous_phase)
            .SetPhase(gas_phase)  // Should throw
            .SetReactants({ CO2 })
            .SetProducts({})
            .SetRateConstant(rate_constant)
            .Build();
      },
      std::system_error);
}

TEST(ChemicalReactionBuilder, SetAerosolScopeAfterSetPhaseThrowsError)
{
  auto CO2 = Species{ "CO2" };
  Phase aqueous_phase{ "aqueous", std::vector<PhaseSpecies>{ CO2 } };
  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ CO2 } };
  auto rate_constant = ArrheniusRateConstant{ { .A_ = 1.0 } };

  EXPECT_THROW(
      {
        ChemicalReactionBuilder()
            .SetPhase(gas_phase)
            .SetAerosolScope("accumulation", aqueous_phase)  // Should throw
            .SetReactants({ CO2 })
            .SetProducts({})
            .SetRateConstant(rate_constant)
            .Build();
      },
      std::system_error);
}

TEST(ChemicalReactionBuilder, UsingOnlySetPhaseWorks)
{
  auto CO2 = Species{ "CO2" };
  auto H2O = Species{ "H2O" };
  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ CO2, H2O } };
  auto rate_constant = ArrheniusRateConstant{ { .A_ = 1.0 } };

  // Should not throw - only SetPhase is used
  Process reaction = ChemicalReactionBuilder()
                         .SetPhase(gas_phase)
                         .SetReactants({ CO2 })
                         .SetProducts({ StoichSpecies(H2O, 1.0) })
                         .SetRateConstant(rate_constant)
                         .Build();

  auto* chem_reaction = std::get_if<ChemicalReaction>(&reaction.process_);
  ASSERT_NE(chem_reaction, nullptr);

  // Verify no scoping was applied
  EXPECT_EQ(chem_reaction->reactants_[0].name_, "CO2");
  EXPECT_EQ(chem_reaction->products_[0].species_.name_, "H2O");
}

TEST(ChemicalReactionBuilder, UsingOnlySetAerosolScopeWorks)
{
  auto CO2 = Species{ "CO2" };
  auto H2O = Species{ "H2O" };
  Phase aqueous_phase{ "aqueous", std::vector<PhaseSpecies>{ CO2, H2O } };
  auto rate_constant = ArrheniusRateConstant{ { .A_ = 1.0 } };

  // Should not throw - only SetAerosolScope is used
  Process reaction = ChemicalReactionBuilder()
                         .SetAerosolScope("accumulation", aqueous_phase)
                         .SetReactants({ CO2 })
                         .SetProducts({ StoichSpecies(H2O, 1.0) })
                         .SetRateConstant(rate_constant)
                         .Build();

  auto* chem_reaction = std::get_if<ChemicalReaction>(&reaction.process_);
  ASSERT_NE(chem_reaction, nullptr);

  // Verify scoping was applied
  EXPECT_EQ(chem_reaction->reactants_[0].name_, "accumulation.aqueous.CO2");
  EXPECT_EQ(chem_reaction->products_[0].species_.name_, "accumulation.aqueous.H2O");
}

// ============================================================================
// Tests for enforcing SetAerosolScope order
// ============================================================================

TEST(ChemicalReactionBuilder, SetAerosolScopeAfterSetReactantsThrowsError)
{
  auto CO2 = Species{ "CO2" };
  Phase aqueous_phase{ "aqueous", std::vector<PhaseSpecies>{ CO2 } };
  auto rate_constant = ArrheniusRateConstant{ { .A_ = 1.0 } };

  EXPECT_THROW(
      {
        ChemicalReactionBuilder()
            .SetReactants({ CO2 })                           // Called first
            .SetAerosolScope("accumulation", aqueous_phase)  // Should throw
            .SetProducts({})
            .SetRateConstant(rate_constant)
            .Build();
      },
      std::system_error);
}

TEST(ChemicalReactionBuilder, SetAerosolScopeAfterSetProductsThrowsError)
{
  auto CO2 = Species{ "CO2" };
  Phase aqueous_phase{ "aqueous", std::vector<PhaseSpecies>{ CO2 } };
  auto rate_constant = ArrheniusRateConstant{ { .A_ = 1.0 } };

  EXPECT_THROW(
      {
        ChemicalReactionBuilder()
            .SetProducts({ StoichSpecies(CO2, 1.0) })                // Called first
            .SetAerosolScope("accumulation", aqueous_phase)  // Should throw
            .SetReactants({})
            .SetRateConstant(rate_constant)
            .Build();
      },
      std::system_error);
}

TEST(ChemicalReactionBuilder, SetAerosolScopeAfterBothReactantsAndProductsThrowsError)
{
  auto CO2 = Species{ "CO2" };
  auto H2O = Species{ "H2O" };
  Phase aqueous_phase{ "aqueous", std::vector<PhaseSpecies>{ CO2, H2O } };
  auto rate_constant = ArrheniusRateConstant{ { .A_ = 1.0 } };

  EXPECT_THROW(
      {
        ChemicalReactionBuilder()
            .SetReactants({ CO2 })
            .SetProducts({ StoichSpecies(H2O, 1.0) })
            .SetAerosolScope("accumulation", aqueous_phase)  // Should throw
            .SetRateConstant(rate_constant)
            .Build();
      },
      std::system_error);
}

TEST(ChemicalReactionBuilder, CorrectOrderSetAerosolScopeBeforeReactantsAndProductsWorks)
{
  auto CO2 = Species{ "CO2" };
  auto H2O = Species{ "H2O" };
  Phase aqueous_phase{ "aqueous", std::vector<PhaseSpecies>{ CO2, H2O } };
  auto rate_constant = ArrheniusRateConstant{ { .A_ = 1.0 } };

  // Should not throw - correct order
  Process reaction = ChemicalReactionBuilder()
                         .SetAerosolScope("accumulation", aqueous_phase)  // Called first
                         .SetReactants({ CO2 })
                         .SetProducts({ StoichSpecies(H2O, 1.0) })
                         .SetRateConstant(rate_constant)
                         .Build();

  auto* chem_reaction = std::get_if<ChemicalReaction>(&reaction.process_);
  ASSERT_NE(chem_reaction, nullptr);

  // Verify scoping was applied correctly
  EXPECT_EQ(chem_reaction->reactants_[0].name_, "accumulation.aqueous.CO2");
  EXPECT_EQ(chem_reaction->products_[0].species_.name_, "accumulation.aqueous.H2O");
}
