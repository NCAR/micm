#include <micm/process/arrhenius_rate_constant.hpp>

#include <gtest/gtest.h>

TEST(ArrheniusRateConstant, DefaultConstructor){
  micm::ArrheniusRateConstant<double> arrhenius{};
}