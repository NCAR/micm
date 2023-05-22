#include <micm/process/photolysis_rate_constant.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

TEST(PhotolysisRateConstant, CalculateWithSystem){
  micm::State state {0,1,1};
  state.custom_rate_parameters_[0] = 0.5;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::PhotolysisRateConstant photo{};
  auto k = photo.calculate(state, params);
  EXPECT_EQ(k, 0.5);
}

TEST(PhotolysisRateConstant, ConstructorWithRate){
  micm::State state {0,1,1};
  state.custom_rate_parameters_[0] = 1.1;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::PhotolysisRateConstant photo{};
  auto k = photo.calculate(state, params);
  EXPECT_EQ(k, 1.1);
}

TEST(PhotolysisRateConstant, ConstructorWithRateAndName){
  micm::State state {0,1,1};
  state.custom_rate_parameters_[0] = 1.1;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::PhotolysisRateConstant photo("a name");
  auto k = photo.calculate(state, params);
  EXPECT_EQ(k, 1.1);
  EXPECT_EQ(photo.name_, "a name");
}