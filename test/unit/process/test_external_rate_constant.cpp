#include <micm/process/external_rate_constant.hpp>

#include <gtest/gtest.h>

TEST(ExternalRateConstant, DefaultConstructor){
  micm::ExternalRateConstant<double> external{};
}