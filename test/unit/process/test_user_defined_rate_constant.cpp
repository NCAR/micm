#include <gtest/gtest.h>

#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

TEST(UserDefinedRateConstant, CalculateWithSystem)
{
  micm::State<micm::Matrix> state{ 0, 1, 1 };
  state.custom_rate_parameters_[0][0] = 0.5;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::UserDefinedRateConstant photo{};
  auto k = photo.calculate(state.conditions_[0], params);
  EXPECT_EQ(k, 0.5);
}

TEST(UserDefinedRateConstant, ConstructorWithRate)
{
  micm::State<micm::Matrix> state{ 0, 1, 1 };
  state.custom_rate_parameters_[0][0] = 1.1;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::UserDefinedRateConstant photo{};
  auto k = photo.calculate(state.conditions_[0], params);
  EXPECT_EQ(k, 1.1);
}

TEST(UserDefinedRateConstant, ConstructorWithRateAndName)
{
  micm::State<micm::Matrix> state{ 0, 1, 1 };
  state.custom_rate_parameters_[0][0] = 1.1;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::UserDefinedRateConstant photo({ .label_ = "a name" });
  auto k = photo.calculate(state.conditions_[0], params);
  EXPECT_EQ(k, 1.1);
  EXPECT_EQ(photo.CustomParameters()[0], "a name");
}
