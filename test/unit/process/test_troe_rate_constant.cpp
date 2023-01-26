#include <micm/process/troe_rate_constant.hpp>

#include <gtest/gtest.h>

TEST(TroeRateConstant, DefaultConstructor){
  micm::TroeRateConstant troe{};
}

TEST(TroeRateConstant, CalculateWithSystem){
  micm::TroeRateConstant troe{};
  auto k = troe.calculate(micm::System());
  ASSERT_TRUE(std::isnan(k));
}

TEST(TroeRateConstant, CalculateWithPrescribedArugments){
  micm::TroeRateConstant troe{};
  auto k = troe.calculate(1.0, 1.0);
  ASSERT_TRUE(std::isnan(k));
}