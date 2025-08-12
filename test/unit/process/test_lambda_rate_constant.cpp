#include <micm/process/rate_constant/lambda_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

TEST(LambdaRateConstant, CalculateWithSystem)
{
  micm::LambdaRateConstant my_lambda{ micm::LambdaRateConstantParameters{
      .label_ = "test lambda",
      .lambda_function_ = [](const micm::Conditions& conditions) { return 2.0 * conditions.temperature_; } } };
  micm::Conditions conditions = {
    .temperature_ = 300.0  // [K]
  };
  auto k = my_lambda.Calculate(conditions);
  EXPECT_NEAR(k, 600.0, 0.00001);

  std::vector<double> custom_parameters = { 1.0, 2.0 };
  auto custom_params_iter = custom_parameters.begin();
  k = my_lambda.Calculate(conditions, custom_params_iter);
  EXPECT_NEAR(k, 600.0, 0.00001);
  EXPECT_EQ(my_lambda.CustomParameters().size(), 0);
}