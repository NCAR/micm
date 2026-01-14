#include <micm/process/transfer_coefficient/henrys_law_coefficient.hpp>
#include <micm/system/conditions.hpp>

#include <gtest/gtest.h>

#include <cmath>

constexpr double TOLERANCE = 1e-13;

TEST(HenrysLawCoefficient, DefaultConstructor)
{
  micm::HenrysLawCoefficient default_hlc{};
  micm::Conditions conditions = {
    .temperature_ = 298.15  // [K]
  };
  conditions.pH = 7.0;

  auto k = default_hlc.Calculate(conditions);
  double expected = 1.0;  // A_ = 1.0, C_ = 0.0, K_a1_ = 0.0, K_a2_ = 0.0
  EXPECT_NEAR(k, expected, TOLERANCE);
}

TEST(HenrysLawCoefficient, CalculateEffective_NeutralPH)
{
  micm::HenrysLawCoefficientParameters parameters;
  parameters.A_ = 32.4;
  parameters.C_ = -2400.0;
  parameters.K_a1_ = 1.0e-7;
  parameters.K_a2_ = 1.0e-11;

  micm::HenrysLawCoefficient hlc(parameters);

  double temperature = 298.15;
  double pH = 7.0;

  double k_h_eff = hlc.CalculateEffective(temperature, pH);
  double k_h = parameters.A_ * std::exp(parameters.C_ / temperature);
  double H_plus = std::pow(10.0, -pH);

  // K_H_eff = K_H * (1 + K_a1/[H+] + K_a1*K_a2/[H+]^2)
  double expected = k_h * (1.0 + parameters.K_a1_ / H_plus + (parameters.K_a1_ * parameters.K_a2_) / (H_plus * H_plus));

  EXPECT_NEAR(k_h_eff, expected, TOLERANCE * expected);
}

TEST(HenrysLawCoefficient, CalculateWithConditions_WithPH)
{
  // Test calculation using Conditions struct with pH
  micm::HenrysLawCoefficientParameters parameters;
  parameters.A_ = 32.4;
  parameters.C_ = -2400.0;
  parameters.K_a1_ = 4.3e-7;
  parameters.K_a2_ = 4.7e-11;

  micm::HenrysLawCoefficient hlc(parameters);

  micm::Conditions conditions = { .temperature_ = 298.15 };
  conditions.pH = 8.0;

  double k = hlc.Calculate(conditions);
  double expected = hlc.CalculateEffective(298.15, 8.0);

  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(HenrysLawCoefficient, CalculateWithConditions_WithoutPH)
{
  // Test calculation using Conditions struct without pH (should default to pH 7.0)
  micm::HenrysLawCoefficientParameters parameters;
  parameters.A_ = 32.4;
  parameters.C_ = -2400.0;
  parameters.K_a1_ = 4.3e-7;
  parameters.K_a2_ = 4.7e-11;

  micm::HenrysLawCoefficient hlc(parameters);

  micm::Conditions conditions = {
    .temperature_ = 298.15  // [K]
  };
  // pH not set

  double k = hlc.Calculate(conditions);
  double expected = hlc.CalculateEffective(298.15, 7.0);  // Default pH = 7.0

  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}

TEST(HenrysLawCoefficient, Clone)
{
  // Test deep copy functionality
  micm::HenrysLawCoefficientParameters parameters;
  parameters.A_ = 32.4;
  parameters.C_ = -2400.0;
  parameters.K_a1_ = 1.0e-7;
  parameters.K_a2_ = 1.0e-11;

  micm::HenrysLawCoefficient original(parameters);
  auto cloned = original.Clone();

  micm::Conditions conditions = { .temperature_ = 298.15 };
  conditions.pH = 7.5;

  double original_result = original.Calculate(conditions);
  double cloned_result = cloned->Calculate(conditions);

  EXPECT_NEAR(original_result, cloned_result, TOLERANCE * original_result);
}

TEST(HenrysLawCoefficient, Calculate_DefaultMethod)
{
  // Test the default Calculate() method (no arguments)
  micm::HenrysLawCoefficientParameters parameters;
  parameters.A_ = 32.4;
  parameters.C_ = -2400.0;
  parameters.K_a1_ = 1.0e-7;
  parameters.K_a2_ = 1.0e-11;

  micm::HenrysLawCoefficient hlc(parameters);

  double k = hlc.Calculate();
  // Should use default T=298.15 K and pH=7.0
  double expected = hlc.CalculateEffective(298.15, 7.0);

  EXPECT_NEAR(k, expected, TOLERANCE * expected);
}
