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
