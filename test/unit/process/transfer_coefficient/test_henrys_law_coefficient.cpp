#include <micm/process/transfer_coefficient/henrys_law_coefficient.hpp>
#include <micm/system/conditions.hpp>
#include <micm/util/constants.hpp>

#include <gtest/gtest.h>

#include <cmath>

constexpr double TOLERANCE = 1e-10;

TEST(HenrysLawConstant, DefaultConstructor)
{
  micm::HenrysLawConstant default_hlc{};
  micm::Conditions conditions = {
    .temperature_ = 298.15  // [K]
  };

  auto k = default_hlc.Calculate(conditions);
  // Default parameters: H_ref_ = 1.3e-3, enthalpy_ = -12000.0, temperature_ref_ = 298.15
  // At T = T_ref, the exponential term equals 1.0, so H = H_ref
  double expected = 1.3e-3;
  EXPECT_NEAR(k, expected, TOLERANCE);
}

TEST(HenrysLawConstant, CalculateAtReferenceTemperature)
{
  micm::HenrysLawConstantParameters parameters;
  parameters.H_ref_ = 1.5e-3;
  parameters.enthalpy_ = -10000.0;
  parameters.temperature_ref_ = 298.15;

  micm::HenrysLawConstant hlc(parameters);

  micm::Conditions conditions = { .temperature_ = 298.15 };

  double k = hlc.Calculate(conditions);
  // At reference temperature, H = H_ref
  double expected = 1.5e-3;

  EXPECT_NEAR(k, expected, TOLERANCE);
}

TEST(HenrysLawConstant, CalculateAtDifferentTemperature)
{
  micm::HenrysLawConstantParameters parameters;
  parameters.H_ref_ = 1.3e-3;
  parameters.enthalpy_ = -12000.0;
  parameters.temperature_ref_ = 298.15;

  micm::HenrysLawConstant hlc(parameters);

  micm::Conditions conditions = { .temperature_ = 288.15 };

  double k = hlc.Calculate(conditions);
  // H = H_ref * exp((-enthalpy / R) * (1/T - 1/T_ref))
  double expected = parameters.H_ref_ * 
                    std::exp((-parameters.enthalpy_ / micm::constants::GAS_CONSTANT) * 
                            (1.0 / 288.15 - 1.0 / parameters.temperature_ref_));

  EXPECT_NEAR(k, expected, TOLERANCE);
}

TEST(HenrysLawConstant, CalculateWithCustomParameters)
{
  // Test with different parameter values
  micm::HenrysLawConstantParameters parameters;
  parameters.H_ref_ = 2.5e-2;
  parameters.enthalpy_ = -15000.0;
  parameters.temperature_ref_ = 300.0;

  micm::HenrysLawConstant hlc(parameters);

  micm::Conditions conditions = { .temperature_ = 310.0 };

  double k = hlc.Calculate(conditions);
  double expected = parameters.H_ref_ * 
                    std::exp((-parameters.enthalpy_ / micm::constants::GAS_CONSTANT) * 
                            (1.0 / 310.0 - 1.0 / 300.0));

  EXPECT_NEAR(k, expected, TOLERANCE * std::abs(expected));
}

TEST(HenrysLawConstant, Clone)
{
  // Test deep copy functionality
  micm::HenrysLawConstantParameters parameters;
  parameters.H_ref_ = 1.8e-3;
  parameters.enthalpy_ = -13500.0;
  parameters.temperature_ref_ = 298.15;

  micm::HenrysLawConstant original(parameters);
  auto cloned = original.Clone();

  micm::Conditions conditions = { .temperature_ = 290.0 };

  double original_result = original.Calculate(conditions);
  double cloned_result = cloned->Calculate(conditions);

  EXPECT_NEAR(original_result, cloned_result, TOLERANCE * std::abs(original_result));
}

TEST(HenrysLawConstant, TemperatureDependence)
{
  // Test that Henry's constant changes correctly with temperature
  micm::HenrysLawConstantParameters parameters;
  parameters.H_ref_ = 1.3e-3;
  parameters.enthalpy_ = -12000.0;  // Negative enthalpy means solubility increases with temperature
  parameters.temperature_ref_ = 298.15;

  micm::HenrysLawConstant hlc(parameters);

  micm::Conditions cold_conditions = { .temperature_ = 280.0 };
  micm::Conditions warm_conditions = { .temperature_ = 320.0 };

  double k_cold = hlc.Calculate(cold_conditions);
  double k_warm = hlc.Calculate(warm_conditions);

  // With negative enthalpy, Henry's constant should increase with temperature
  EXPECT_GT(k_warm, k_cold);
}
