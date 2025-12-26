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

TEST(HenrysLawCoefficient, CalculateKH_SimpleCase)
{
  // Test basic Henry's Law constant calculation: K_H = A * exp(C/T)
  micm::HenrysLawCoefficientParameters parameters;
  parameters.A_ = 32.4;
  parameters.C_ = -3.2e-3;

  micm::HenrysLawCoefficient hlc(parameters);
  
  double temperature = 298.15;  // [K]
  double k_h = hlc.CalculateKH(temperature);
  double expected = 32.4 * std::exp(-3.2e-3 / 298.15);
  EXPECT_NEAR(k_h, expected, TOLERANCE * expected);
}

TEST(HenrysLawCoefficient, CalculateKH_TemperatureDependence)
{
  // Test temperature dependence
  micm::HenrysLawCoefficientParameters parameters;
  parameters.A_ = 32.4;
  parameters.C_ = -2400.0;  // More typical value for C parameter [K]

  micm::HenrysLawCoefficient hlc(parameters);
  
  // Test at different temperatures
  double T1 = 273.15;  // 0°C
  double T2 = 298.15;  // 25°C
  double T3 = 313.15;  // 40°C

  double k_h1 = hlc.CalculateKH(T1);
  double k_h2 = hlc.CalculateKH(T2);
  double k_h3 = hlc.CalculateKH(T3);

  double expected1 = 32.4 * std::exp(-2400.0 / T1);
  double expected2 = 32.4 * std::exp(-2400.0 / T2);
  double expected3 = 32.4 * std::exp(-2400.0 / T3);

  EXPECT_NEAR(k_h1, expected1, TOLERANCE * expected1);
  EXPECT_NEAR(k_h2, expected2, TOLERANCE * expected2);
  EXPECT_NEAR(k_h3, expected3, TOLERANCE * expected3);

  // Verify that Henry's Law constant decreases with temperature for negative C
  EXPECT_GT(k_h1, k_h2);
  EXPECT_GT(k_h2, k_h3);
}

TEST(HenrysLawCoefficient, CalculateEffective_NeutralPH)
{
  // Test effective Henry's Law coefficient at neutral pH
  micm::HenrysLawCoefficientParameters parameters;
  parameters.A_ = 32.4;
  parameters.C_ = -2400.0;
  parameters.K_a1_ = 1.0e-7;  // pK_a1 ~ 7
  parameters.K_a2_ = 1.0e-11; // pK_a2 ~ 11

  micm::HenrysLawCoefficient hlc(parameters);
  
  double temperature = 298.15;  // [K]
  double pH = 7.0;
  
  double k_h_eff = hlc.CalculateEffective(temperature, pH);
  double k_h = hlc.CalculateKH(temperature);
  double H_plus = std::pow(10.0, -pH);
  
  // K_H_eff = K_H * (1 + K_a1/[H+] + K_a1*K_a2/[H+]^2)
  double expected = k_h * (1.0 + parameters.K_a1_ / H_plus + 
                           (parameters.K_a1_ * parameters.K_a2_) / (H_plus * H_plus));
  
  EXPECT_NEAR(k_h_eff, expected, TOLERANCE * expected);
}

TEST(HenrysLawCoefficient, CalculateEffective_AcidicPH)
{
  // Test at acidic pH where dissociation is suppressed
  micm::HenrysLawCoefficientParameters parameters;
  parameters.A_ = 32.4;
  parameters.C_ = -2400.0;
  parameters.K_a1_ = 4.3e-7;   // CO2: pK_a1 ~ 6.37
  parameters.K_a2_ = 4.7e-11;  // CO2: pK_a2 ~ 10.33

  micm::HenrysLawCoefficient hlc(parameters);
  
  double temperature = 298.15;  // [K]
  double pH = 3.0;  // Acidic
  
  double k_h_eff = hlc.CalculateEffective(temperature, pH);
  double k_h = hlc.CalculateKH(temperature);
  
  // At low pH, K_H_eff should be close to K_H (minimal dissociation)
  EXPECT_NEAR(k_h_eff / k_h, 1.0, 0.01);  // Within 1%
}

TEST(HenrysLawCoefficient, CalculateEffective_BasicPH)
{
  // Test at basic pH where dissociation is significant
  micm::HenrysLawCoefficientParameters parameters;
  parameters.A_ = 32.4;
  parameters.C_ = -2400.0;
  parameters.K_a1_ = 4.3e-7;   // CO2: pK_a1 ~ 6.37
  parameters.K_a2_ = 4.7e-11;  // CO2: pK_a2 ~ 10.33

  micm::HenrysLawCoefficient hlc(parameters);
  
  double temperature = 298.15;  // [K]
  double pH = 10.0;  // Basic
  
  double k_h_eff = hlc.CalculateEffective(temperature, pH);
  double k_h = hlc.CalculateKH(temperature);
  double H_plus = std::pow(10.0, -pH);
  
  // K_H_eff = K_H * (1 + K_a1/[H+] + K_a1*K_a2/[H+]^2)
  double expected = k_h * (1.0 + parameters.K_a1_ / H_plus + 
                           (parameters.K_a1_ * parameters.K_a2_) / (H_plus * H_plus));
  
  EXPECT_NEAR(k_h_eff, expected, TOLERANCE * expected);
  
  // At high pH, K_H_eff should be much larger than K_H
  EXPECT_GT(k_h_eff, k_h * 10.0);
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
  
  micm::Conditions conditions = {
    .temperature_ = 298.15  // [K]
  };
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

TEST(HenrysLawCoefficient, CO2_RealisticParameters)
{
  // Test with realistic CO2 parameters
  // Based on comment: CO2(g) <-> H2CO3(aq) <-> HCO3- + H+ <-> CO32- + 2H+
  micm::HenrysLawCoefficientParameters parameters;
  parameters.A_ = 1.14e-2;
  parameters.C_ = 2300.0;
  parameters.K_a1_ = 1.0e-5;
  parameters.K_a2_ = 2.0e-5;

  micm::HenrysLawCoefficient hlc(parameters);
  
  // Test at seawater conditions (T ~ 288 K, pH ~ 8.1)
  double temperature = 288.15;
  double pH = 8.1;
  
  double k_h = hlc.CalculateKH(temperature);
  double expected_k_h = 1.14e-2 * std::exp(2300.0 / temperature);
  EXPECT_NEAR(k_h, expected_k_h, TOLERANCE * expected_k_h);
  
  double k_h_eff = hlc.CalculateEffective(temperature, pH);
  double H_plus = std::pow(10.0, -pH);
  double expected_k_h_eff = k_h * (1.0 + parameters.K_a1_ / H_plus + 
                                   (parameters.K_a1_ * parameters.K_a2_) / (H_plus * H_plus));
  EXPECT_NEAR(k_h_eff, expected_k_h_eff, TOLERANCE * expected_k_h_eff);
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
  
  micm::Conditions conditions = {
    .temperature_ = 298.15
  };
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
