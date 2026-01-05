#include <micm/process/rate_constant/reversible_rate_constant.hpp>
#include <micm/system/conditions.hpp>

#include <gtest/gtest.h>

#include <cmath>

constexpr double TOLERANCE = 1e-13;

TEST(ReversibleRateConstant, DefaultConstructor)
{
  micm::ReversibleRateConstant rc{};
  micm::Conditions conditions = { .temperature_ = 298.15 };

  auto k = rc.Calculate(conditions);
  // Default: A_ = 1, C_ = 0, k_r_ = 0
  // k_f = K_eq * k_r = (1 * exp(0/T)) * 0 = 0
  double expected = 0.0;
  EXPECT_NEAR(k, expected, TOLERANCE);
}

TEST(ReversibleRateConstant, CalculateWithSystem)
{
  micm::Conditions conditions = { .temperature_ = 301.24 };

  // Test with specific parameters
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.14e-2;
  parameters.C_ = 2300.0;
  parameters.k_r_ = 0.32;

  micm::ReversibleRateConstant rc(parameters);
  auto k = rc.Calculate(conditions);

  // k_f = K_eq * k_r = (A * exp(C/T)) * k_r
  double K_eq = 1.14e-2 * std::exp(2300.0 / 301.24);
  double expected = K_eq * 0.32;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(ReversibleRateConstant, CalculateWithPrescribedArguments)
{
  double temperature = 298.15;

  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.14e-2;
  parameters.C_ = 2300.0;
  parameters.k_r_ = 0.32;

  micm::ReversibleRateConstant rc(parameters);
  auto k = rc.Calculate(temperature);

  // k_f = K_eq * k_r = (A * exp(C/T)) * k_r
  double K_eq = 1.14e-2 * std::exp(2300.0 / temperature);
  double expected = K_eq * 0.32;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(ReversibleRateConstant, TemperatureDependence)
{
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.5;
  parameters.C_ = 2000.0;
  parameters.k_r_ = 0.5;

  micm::ReversibleRateConstant rc(parameters);

  double T1 = 273.15;  // Low temperature
  double T2 = 298.15;  // Room temperature
  double T3 = 323.15;  // High temperature

  double k1 = rc.Calculate(T1);
  double k2 = rc.Calculate(T2);
  double k3 = rc.Calculate(T3);

  double K_eq1 = 1.5 * std::exp(2000.0 / T1);
  double K_eq2 = 1.5 * std::exp(2000.0 / T2);
  double K_eq3 = 1.5 * std::exp(2000.0 / T3);

  double expected1 = K_eq1 * 0.5;
  double expected2 = K_eq2 * 0.5;
  double expected3 = K_eq3 * 0.5;

  EXPECT_NEAR(k1, expected1, TOLERANCE * expected1);
  EXPECT_NEAR(k2, expected2, TOLERANCE * expected2);
  EXPECT_NEAR(k3, expected3, TOLERANCE * expected3);
}

TEST(ReversibleRateConstant, ZeroReverseRate)
{
  // Test edge case where k_r = 0 (effectively irreversible)
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 5.0;
  parameters.C_ = 1000.0;
  parameters.k_r_ = 0.0;

  micm::ReversibleRateConstant rc(parameters);

  micm::Conditions conditions = { .temperature_ = 298.15 };

  auto k = rc.Calculate(conditions);
  EXPECT_NEAR(k, 0.0, TOLERANCE);
}

TEST(ReversibleRateConstant, EquilibriumConstantOnly)
{
  // Test case where k_r = 1.0, so k_f = K_eq
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 3.5;
  parameters.C_ = 1800.0;
  parameters.k_r_ = 1.0;

  micm::ReversibleRateConstant rc(parameters);

  double temperature = 298.15;
  auto k = rc.Calculate(temperature);

  double K_eq = 3.5 * std::exp(1800.0 / temperature);
  EXPECT_NEAR(k, K_eq, TOLERANCE * K_eq);
}

TEST(ReversibleRateConstant, Clone)
{
  // Test deep copy functionality
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.14e-2;
  parameters.C_ = 2300.0;
  parameters.k_r_ = 0.32;

  micm::ReversibleRateConstant original(parameters);
  auto cloned = original.Clone();

  micm::Conditions conditions = { .temperature_ = 298.15 };

  double original_result = original.Calculate(conditions);
  double cloned_result = cloned->Calculate(conditions);

  EXPECT_NEAR(original_result, cloned_result, TOLERANCE * original_result);
}