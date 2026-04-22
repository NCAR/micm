// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/process/reaction_rate_store.hpp>
#include <micm/system/system.hpp>
#include <micm/util/constants.hpp>

#include <gtest/gtest.h>

#define _USE_MATH_DEFINES
#include <math.h>
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

using namespace micm;

constexpr double TOLERANCE = 1e-13;

TEST(SurfaceRateConstant, CalculateDefaultProbability)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  double foo_diffusion_coefficient = 2.3e2;
  PhaseSpecies foo_gas_species(foo, foo_diffusion_coefficient);

  // Build SurfaceRateConstantData as BuildFrom would
  SurfaceRateConstantData data;
  data.diffusion_coefficient_  = foo_diffusion_coefficient;
  data.mean_free_speed_factor_ = 8.0 * constants::GAS_CONSTANT / (M_PI * 0.025);
  data.reaction_probability_   = 1.0;  // default
  data.custom_param_base_index_ = 0;

  double temperature = 273.65;  // K
  double custom_params[2] = { 1.0e-7, 2.5e6 };  // effective radius [m], particle number concentration [# m-3]
  double k;
  CalculateSurface(&data, 1, temperature, custom_params, &k);

  double expected = 4.0 * 2.5e6 * M_PI * std::pow(1.0e-7, 2) /
                    (1.0e-7 / 2.3e2 + 4.0 / (std::sqrt(8.0 * constants::GAS_CONSTANT * 273.65 / (M_PI * 0.025))));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(SurfaceRateConstant, CalculateSpecifiedProbability)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  double foo_diffusion_coefficient = 2.3e2;
  PhaseSpecies foo_gas_species(foo, foo_diffusion_coefficient);

  SurfaceRateConstantData data;
  data.diffusion_coefficient_  = foo_diffusion_coefficient;
  data.mean_free_speed_factor_ = 8.0 * constants::GAS_CONSTANT / (M_PI * 0.025);
  data.reaction_probability_   = 0.74;
  data.custom_param_base_index_ = 0;

  double temperature = 273.65;
  double custom_params[2] = { 1.0e-7, 2.5e6 };
  double k;
  CalculateSurface(&data, 1, temperature, custom_params, &k);

  double expected = 4.0 * 2.5e6 * M_PI * std::pow(1.0e-7, 2) /
                    (1.0e-7 / 2.3e2 + 4.0 / (std::sqrt(8.0 * constants::GAS_CONSTANT * 273.65 / (M_PI * 0.025)) * 0.74));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(SurfaceRateConstant, DiffusionCoefficientIsMissing)
{
  // BuildFrom throws when diffusion_coefficient is missing
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 } });
  PhaseSpecies foo_gas_species(foo);  // no diffusion coefficient

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ foo_gas_species } };
  SurfaceRateConstantParameters params{ .label_ = "foo", .phase_species_ = foo_gas_species };
  Process proc = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(params).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { proc };

  try
  {
    auto store = ReactionRateConstantStore::BuildFrom(processes);
    FAIL() << "Expected MicmException to be thrown";
  }
  catch (const micm::MicmException& e)
  {
    EXPECT_EQ(e.code_, MICM_SPECIES_ERROR_CODE_PROPERTY_NOT_FOUND);
    EXPECT_NE(std::string(e.what()).find("Diffusion coefficient for species 'foo' is not defined"), std::string::npos);
  }
}

TEST(SurfaceRateConstant, MolecularWeightIsMissing)
{
  // BuildFrom throws when molecular weight is missing (GetProperty throws for missing property)
  Species foo("foo");
  double foo_diffusion_coefficient = 2.3e2;
  PhaseSpecies foo_gas_species(foo, foo_diffusion_coefficient);

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ foo_gas_species } };
  SurfaceRateConstantParameters params{ .label_ = "foo", .phase_species_ = foo_gas_species };
  Process proc = ChemicalReactionBuilder().SetReactants({ foo }).SetRateConstant(params).SetPhase(gas_phase).Build();
  std::vector<Process> processes = { proc };

  try
  {
    auto store = ReactionRateConstantStore::BuildFrom(processes);
    FAIL() << "Expected MicmException to be thrown";
  }
  catch (const micm::MicmException& e)
  {
    EXPECT_EQ(e.code_, MICM_SPECIES_ERROR_CODE_PROPERTY_NOT_FOUND);
    // GetProperty throws: "Species: 'foo' Property: 'molecular weight [kg mol-1]'"
    EXPECT_NE(std::string(e.what()).find("foo"), std::string::npos);
    EXPECT_NE(std::string(e.what()).find("molecular weight"), std::string::npos);
  }
}
