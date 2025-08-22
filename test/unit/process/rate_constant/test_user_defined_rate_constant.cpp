#include <micm/process/rate_constant/user_defined_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

TEST(UserDefinedRateConstant, CalculateWithSystem)
{
  auto state_parameters_ = micm::StateParameters{
    .number_of_rate_constants_ = 1,
    .variable_names_ = {"user"},
    .custom_rate_parameter_labels_ = { "my rate", },
  };

  micm::State state{ state_parameters_, 1 };
  state.custom_rate_parameters_[0][0] = 0.5;

  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::UserDefinedRateConstant photo{};
  auto k = photo.Calculate(state.conditions_[0], params);
  EXPECT_EQ(k, 0.5);
}

TEST(UserDefinedRateConstant, ConstructorWithRate)
{
  auto state_parameters_ = micm::StateParameters{
    .number_of_rate_constants_ = 1,
    .variable_names_ = {"user"},
    .custom_rate_parameter_labels_ = { "my rate", },
  };

  micm::State state{ state_parameters_, 1 };
  state.custom_rate_parameters_[0][0] = 1.1;

  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::UserDefinedRateConstant photo{};
  auto k = photo.Calculate(state.conditions_[0], params);
  EXPECT_EQ(k, 1.1);
}

TEST(UserDefinedRateConstant, ConstructorWithRateAndName)
{
  auto state_parameters_ = micm::StateParameters{
    .number_of_rate_constants_ = 1,
    .variable_names_ = {"user"},
    .custom_rate_parameter_labels_ = { "my rate", },
  };

  micm::State state{ state_parameters_, 1 };
  state.custom_rate_parameters_[0][0] = 1.1;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::UserDefinedRateConstant photo({ .label_ = "a name" });
  auto k = photo.Calculate(state.conditions_[0], params);
  EXPECT_EQ(k, 1.1);
  EXPECT_EQ(photo.CustomParameters()[0], "a name");
}

TEST(UserDefinedRateConstant, CustomScalingFactor)
{
  auto state_parameters = micm::StateParameters{
    .number_of_rate_constants_ = 1,
    .variable_names_ = {"user"},
    .custom_rate_parameter_labels_ = { "my rate", },
  };

  micm::State state{ state_parameters, 1 };
  state.custom_rate_parameters_[0][0] = 1.2;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::UserDefinedRateConstant photo({ .label_ = "a name", .scaling_factor_ = 2.0 });
  auto k = photo.Calculate(state.conditions_[0], params);
  EXPECT_EQ(k, 1.2 * 2.0);
  EXPECT_EQ(photo.CustomParameters()[0], "a name");
}
