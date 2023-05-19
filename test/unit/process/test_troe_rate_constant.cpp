#include <micm/process/troe_rate_constant.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

//
// TODO: jiwon 5/18 - comment out because TroeRateConstant and System dosn't have default constructor anymore
//
// TEST(TroeRateConstant, DefaultConstructor){
//   micm::TroeRateConstant troe{};
// }

// TEST(TroeRateConstant, CalculateWithSystem){
//   micm::TroeRateConstant troe{};
//   auto k = troe.calculate(micm::System());
//   ASSERT_TRUE(std::isnan(k));
// }

// TEST(TroeRateConstant, CalculateWithPrescribedArugments){
//   micm::TroeRateConstant troe{};
//   auto k = troe.calculate(1.0, 1.0);
//   ASSERT_TRUE(std::isnan(k));
// }