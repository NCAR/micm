#include <micm/process/troe_rate_constant.hpp>

#include <gtest/gtest.h>
#include <cassert>

TEST(TroeRateConstant, DefaultConstructor){
  micm::TroeRateConstant troe{};
}

TEST(TroeRateConstant, CalculateWithSystem){
  micm::TroeRateConstant troe{};
  auto k = troe.calculate(micm::System());
  assert(std::isnan(k));
}

TEST(TroeRateConstant, CalculateWithPrescribedArugments){
  micm::TroeRateConstant troe{};
  auto k = troe.calculate(1.0, 1.0);
  assert(std::isnan(k));
}