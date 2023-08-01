#include <gtest/gtest.h>

#include <micm/process/surface_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

TEST(SurfaceRateConstant, CalculateDefaultProbability)
{
  micm::Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 }, { "diffusion coefficient [m2 s-1]", 2.3e2 } });
  micm::State<micm::Matrix> state{ 0, 2, 1 };
  state.custom_rate_parameters_[0][0] = 1.0e-7;  // effective radius [m-1]
  state.custom_rate_parameters_[0][1] = 2.5e6;   // particle concentration [# m-3]
  state.conditions_[0].temperature_ = 273.65;    // K
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::SurfaceRateConstant surface{ { .label_ = "foo", .species_ = foo } };
  auto k = surface.calculate(state.conditions_[0], params);
  double k_test = 4.0 * 2.5e6 * M_PI * std::pow(1.0e-7, 2) /
                  (1.0e-7 / 2.3e2 + 4.0 / (std::sqrt(8.0 * GAS_CONSTANT * 273.65 / (M_PI * 0.025))));
  EXPECT_NEAR(k, k_test, k_test * 1.0e-10);
}

TEST(SurfaceRateConstant, CalculateSpecifiedProbability)
{
  micm::Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 }, { "diffusion coefficient [m2 s-1]", 2.3e2 } });
  micm::State<micm::Matrix> state{ 0, 2, 1 };
  state.custom_rate_parameters_[0][0] = 1.0e-7;  // effective radius [m-1]
  state.custom_rate_parameters_[0][1] = 2.5e6;   // particle concentration [# m-3]
  state.conditions_[0].temperature_ = 273.65;    // K
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::SurfaceRateConstant surface{ { .label_ = "foo", .species_ = foo, .reaction_probability_ = 0.74 } };
  auto k = surface.calculate(state.conditions_[0], params);
  double k_test = 4.0 * 2.5e6 * M_PI * std::pow(1.0e-7, 2) /
                  (1.0e-7 / 2.3e2 + 4.0 / (std::sqrt(8.0 * GAS_CONSTANT * 273.65 / (M_PI * 0.025)) * 0.74));
  EXPECT_NEAR(k, k_test, k_test * 1.0e-10);
}