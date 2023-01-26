#include <micm/process/arrhenius_rate_constant.hpp>

#include <gtest/gtest.h>
#include <cassert>

TEST(ArrheniusRateConstant, DefaultConstructor){
  micm::ArrheniusRateConstant arrhenius{};
}

TEST(ArrheniusRateConstant, CalculateWithSystem){
  micm::ArrheniusRateConstant zero{};
  auto k = zero.calculate(micm::System());
  assert(k == 0);

  micm::ArrheniusRateConstant basic(1, 0, 0, 0, 0);
  k = basic.calculate(micm::System());
  assert(k == 0);
}

TEST(ArrheniusRateConstant, CalculateWithPrescribedArugments){
  micm::ArrheniusRateConstant zero{};
  auto k = zero.calculate(1.0, 1.0);
  assert(k == 0);

  micm::ArrheniusRateConstant basic(1, 0, 0, 0, 0);
  k = basic.calculate(micm::System());
  assert(k == 0);
}