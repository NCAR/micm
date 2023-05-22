#include <micm/process/troe_rate_constant.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

TEST(TroeRateConstant, CalculateWithPrescribedArugments){
  micm::State state {0,0,1};
  state.temperature_ = 301.24; // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::TroeRateConstant troe{};
  auto k = troe.calculate(state, params);
  ASSERT_TRUE(std::isnan(k));
}