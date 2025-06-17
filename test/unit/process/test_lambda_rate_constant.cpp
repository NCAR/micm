#include <micm/process/lambda_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

TEST(LambdaRateConstant, CalculateWithSystem)
{
  micm::LambdaRateConstant my_lambda{
      micm::LamdaRateConstantParameters{
          .label_ = "test lambda",
          .lambda_function_ = [](const micm::Conditions& conditions) { return 2.0 * conditions.temperature_; }
      }
  };
  micm::Conditions conditions = {
    .temperature_ = 300.0  // [K]
  };
  auto k = my_lambda.Calculate(conditions);
  EXPECT_NEAR(k, 600.0, 0.00001);
}