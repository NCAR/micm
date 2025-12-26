#include <micm/process/rate_constant/reversible_rate_constant.hpp>
#include <micm/system/conditions.hpp>

#include <gtest/gtest.h>

#include <cmath>

constexpr double TOLERANCE = 1e-13;

TEST(ReversibleRateConstant, DefaultConstructor)
{
  micm::ReversibleRateConstant default_rc{};
  micm::Conditions conditions = {
    .temperature_ = 298.15  // [K]
  };

  auto k = default_rc.Calculate(conditions);
  // Default: A_ = 1, C_ = 0, k_r_ = 0
  // k_f = K_eq * k_r = (1 * exp(0/T)) * 0 = 0
  double expected = 0.0;
  EXPECT_NEAR(k, expected, TOLERANCE);
}

TEST(ReversibleRateConstant, CalculateWithSystem)
{
  micm::Conditions conditions = {
    .temperature_ = 301.24  // [K]
  };

  // Test with specific parameters from the API example
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
  double temperature = 298.15;  // [K]
  
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
  // Test that rate constant changes with temperature as expected
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.5;
  parameters.C_ = 2000.0;  // Positive C means K_eq increases with temperature
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
  
  // For positive C, K_eq decreases with temperature, so k_f should also decrease
  EXPECT_GT(k1, k2);
  EXPECT_GT(k2, k3);
}

TEST(ReversibleRateConstant, NegativeActivationEnergy)
{
  // Test with negative C (common for some reactions)
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 2.0;
  parameters.C_ = -1500.0;  // Negative C means K_eq decreases with temperature
  parameters.k_r_ = 0.25;

  micm::ReversibleRateConstant rc(parameters);
  
  double T1 = 273.15;
  double T2 = 323.15;
  
  double k1 = rc.Calculate(T1);
  double k2 = rc.Calculate(T2);
  
  // For negative C, K_eq increases with temperature
  EXPECT_LT(k1, k2);
}

TEST(ReversibleRateConstant, ZeroReverseRate)
{
  // Test edge case where k_r = 0 (effectively irreversible)
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 5.0;
  parameters.C_ = 1000.0;
  parameters.k_r_ = 0.0;

  micm::ReversibleRateConstant rc(parameters);
  
  micm::Conditions conditions = {
    .temperature_ = 298.15
  };
  
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

TEST(ReversibleRateConstant, CalculateWithConditionsAndCustomParameters)
{
  // Test the Calculate method that takes custom_parameters iterator
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.14e-2;
  parameters.C_ = 2300.0;
  parameters.k_r_ = 0.32;

  micm::ReversibleRateConstant rc(parameters);
  
  micm::Conditions conditions = {
    .temperature_ = 301.24
  };
  
  std::vector<double> custom_params = {1.0, 2.0, 3.0};
  auto k = rc.Calculate(conditions, custom_params.begin());
  
  // Should ignore custom parameters and use only temperature
  double K_eq = 1.14e-2 * std::exp(2300.0 / 301.24);
  double expected = K_eq * 0.32;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
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
  
  micm::Conditions conditions = {
    .temperature_ = 298.15
  };
  
  double original_result = original.Calculate(conditions);
  double cloned_result = cloned->Calculate(conditions);
  
  EXPECT_NEAR(original_result, cloned_result, TOLERANCE * original_result);
}

TEST(ReversibleRateConstant, LargeEquilibriumConstant)
{
  // Test with large K_eq (highly favorable forward reaction)
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.0e10;
  parameters.C_ = 5000.0;
  parameters.k_r_ = 0.01;

  micm::ReversibleRateConstant rc(parameters);
  
  double temperature = 300.0;
  auto k = rc.Calculate(temperature);
  
  double K_eq = 1.0e10 * std::exp(5000.0 / temperature);
  double expected = K_eq * 0.01;
  EXPECT_NEAR(k, expected, TOLERANCE * expected);
  
  // Verify the result is physically reasonable (very large)
  EXPECT_GT(k, 1.0e10);
}

TEST(ReversibleRateConstant, SmallEquilibriumConstant)
{
  // Test with small K_eq (unfavorable forward reaction)
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.0e-10;
  parameters.C_ = -5000.0;
  parameters.k_r_ = 0.1;

  micm::ReversibleRateConstant rc(parameters);
  
  double temperature = 300.0;
  auto k = rc.Calculate(temperature);
  
  double K_eq = 1.0e-10 * std::exp(-5000.0 / temperature);
  double expected = K_eq * 0.1;
  EXPECT_NEAR(k, expected, TOLERANCE * std::abs(expected));
  
  // Verify the result is physically reasonable (very small)
  EXPECT_LT(k, 1.0e-20);
}

TEST(ReversibleRateConstant, StandardConditions)
{
  // Test at standard conditions (298.15 K, 1 atm)
  micm::ReversibleRateConstantParameters parameters;
  parameters.A_ = 1.0;
  parameters.C_ = 0.0;
  parameters.k_r_ = 1.0;

  micm::ReversibleRateConstant rc(parameters);
  
  micm::Conditions conditions = {
    .temperature_ = 298.15,
    .pressure_ = 101325.0  // 1 atm in Pa
  };
  
  auto k = rc.Calculate(conditions);
  // K_eq = 1.0 * exp(0) = 1.0
  // k_f = 1.0 * 1.0 = 1.0
  double expected = 1.0;
  EXPECT_NEAR(k, expected, TOLERANCE);
}
