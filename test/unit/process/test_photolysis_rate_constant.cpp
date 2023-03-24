#include <micm/process/photolysis_rate_constant.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

TEST(PhotolysisRateConstant, DefaultConstructor){
  micm::PhotolysisRateConstant rate{};
}

TEST(PhotolysisRateConstant, CalculateWithSystem){
  micm::PhotolysisRateConstant rate{0.5};
  auto k = rate.calculate(micm::System());
  EXPECT_EQ(k, 0.5);
}

TEST(PhotolysisRateConstant, ConstructorWithRate){
  micm::PhotolysisRateConstant rate(1.1);
  auto k = rate.calculate(micm::System());
  EXPECT_EQ(k, 1.1);
}

TEST(PhotolysisRateConstant, ConstructorWithRateAndName){
  micm::PhotolysisRateConstant rate(1.1, "a name");
  auto k = rate.calculate(micm::System());
  EXPECT_EQ(k, 1.1);
  EXPECT_EQ(rate.name_, "a name");
}