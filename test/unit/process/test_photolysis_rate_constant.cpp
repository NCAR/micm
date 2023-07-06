#include <gtest/gtest.h>

#include <micm/process/photolysis_rate_constant.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>

TEST(PhotolysisRateConstant, CalculateWithSystem)
{
  micm::State<micm::Matrix> state{ 0, 1, 1 };
  state.custom_rate_parameters_[0][0] = 0.5;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::PhotolysisRateConstant photo{};
  auto k = photo.calculate(state.conditions_[0], params);
  EXPECT_EQ(k, 0.5);
}

TEST(PhotolysisRateConstant, ConstructorWithRate)
{
  micm::State<micm::Matrix> state{ 0, 1, 1 };
  state.custom_rate_parameters_[0][0] = 1.1;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::PhotolysisRateConstant photo{};
  auto k = photo.calculate(state.conditions_[0], params);
  EXPECT_EQ(k, 1.1);
}

TEST(PhotolysisRateConstant, ConstructorWithRateAndName)
{
  micm::State<micm::Matrix> state{ 0, 1, 1 };
  state.custom_rate_parameters_[0][0] = 1.1;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::PhotolysisRateConstant photo("a name");
  auto k = photo.calculate(state.conditions_[0], params);
  EXPECT_EQ(k, 1.1);
  EXPECT_EQ(photo.name_, "a name");
}
