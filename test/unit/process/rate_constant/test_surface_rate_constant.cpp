#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

#include <memory>

using namespace micm;

constexpr double TOLERANCE = 1e-13;

TEST(SurfaceRateConstant, CalculateDefaultProbability)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 }});
  double foo_diffusion_coefficient = 2.3e2;
  PhaseSpecies foo_gas_species(foo, foo_diffusion_coefficient);

  auto state_parameters_ = StateParameters{
    .number_of_rate_constants_ = 1,
    .variable_names_ = { "surface" },
    .custom_rate_parameter_labels_ = { "effective radius [m]", "particle number concentration [# m-3]" },
  };
  State state{ state_parameters_, 1 };
  state.custom_rate_parameters_[0][0] = 1.0e-7;  // effective radius [m]
  state.custom_rate_parameters_[0][1] = 2.5e6;   // particle concentration [# m-3]
  state.conditions_[0].temperature_ = 273.65;    // K
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();

  SurfaceRateConstantParameters parameters{
    .label_ = "foo",
    .phase_species_ = foo_gas_species
  };
  SurfaceRateConstant surface(parameters);

  EXPECT_EQ(surface.SizeCustomParameters(), 2);
  EXPECT_EQ(surface.CustomParameters()[0], "foo.effective radius [m]");
  EXPECT_EQ(surface.CustomParameters()[1], "foo.particle number concentration [# m-3]");
  auto k = surface.Calculate(state.conditions_[0], params);
  double expected = 4.0 * 2.5e6 * M_PI * std::pow(1.0e-7, 2) /
                    (1.0e-7 / 2.3e2 + 4.0 / (std::sqrt(8.0 * constants::GAS_CONSTANT * 273.65 / (M_PI * 0.025))));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(SurfaceRateConstant, CalculateSpecifiedProbability)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 }});
  double foo_diffusion_coefficient = 2.3e2;
  PhaseSpecies foo_gas_species(foo, foo_diffusion_coefficient);

  auto state_parameters_ = StateParameters{
    .number_of_rate_constants_ = 1,
    .variable_names_ = { "surface" },
    .custom_rate_parameter_labels_ = { "effective radius [m]", "particle number concentration [# m-3]" },
  };
  State state{ state_parameters_, 1 };
  state.custom_rate_parameters_[0][0] = 1.0e-7;  // effective radius [m]
  state.custom_rate_parameters_[0][1] = 2.5e6;   // particle concentration [# m-3]
  state.conditions_[0].temperature_ = 273.65;    // K
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  SurfaceRateConstantParameters parameters{
    .label_ = "foo",
    .phase_species_ = foo_gas_species,
    .reaction_probability_ = 0.74
  };
  SurfaceRateConstant surface(parameters);

  EXPECT_EQ(surface.SizeCustomParameters(), 2);
  EXPECT_EQ(surface.CustomParameters()[0], "foo.effective radius [m]");
  EXPECT_EQ(surface.CustomParameters()[1], "foo.particle number concentration [# m-3]");
  auto k = surface.Calculate(state.conditions_[0], params);
  double expected =
      4.0 * 2.5e6 * M_PI * std::pow(1.0e-7, 2) /
      (1.0e-7 / 2.3e2 + 4.0 / (std::sqrt(8.0 * constants::GAS_CONSTANT * 273.65 / (M_PI * 0.025)) * 0.74));
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(SurfaceRateConstant, DiffusionCoefficientIsMissing)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 }});
  PhaseSpecies foo_gas_species(foo);
  SurfaceRateConstantParameters parameters{
    .label_ = "foo",
    .phase_species_ = foo_gas_species
  };

  try
  {
    SurfaceRateConstant surface_rate_constant(parameters);
  }
  catch (const std::system_error& e)
  {
    EXPECT_EQ(e.code(), make_error_code(MicmSpeciesErrc::PropertyNotFound));
    EXPECT_NE(std::string(e.what()).find("Diffusion coefficient for species 'foo' is not defined"), std::string::npos);
    return;
  }
  FAIL() << "Expected std::system_error to be thrown";
}

TEST(SurfaceRateConstant, MolecularWeightIsMissing)
{
  Species foo("foo");
  double foo_diffusion_coefficient = 2.3e2;
  PhaseSpecies foo_gas_species(foo, foo_diffusion_coefficient);
  SurfaceRateConstantParameters parameters{
    .label_ = "foo",
    .phase_species_ = foo_gas_species
  };

  try
  {
    SurfaceRateConstant surface_rate_constant(parameters);
  }
  catch (const std::system_error& e)
  {
    EXPECT_EQ(e.code(), make_error_code(MicmSpeciesErrc::PropertyNotFound));
    EXPECT_NE(std::string(e.what()).find("Molecular weight for species 'foo' is not defined"), std::string::npos);
    return;
  }
  FAIL() << "Expected std::system_error to be thrown";
}