#include <micm/process/intraphase_process.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/external_rate_constant.hpp>
#include <micm/system/species.hpp>

#include <gtest/gtest.h>

TEST(IntraPhaseProcess, DefaultConstructor){
  micm::IntraPhaseProcess<micm::ArrheniusRateConstant> arrhenius;
  micm::IntraPhaseProcess<micm::TroeRateConstant> troe;
  micm::IntraPhaseProcess<micm::ExternalRateConstant> external;
}