#include <micm/process/intraphase_process.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/external_rate_constant.hpp>
#include <micm/system/species.hpp>

#include <gtest/gtest.h>

TEST(IntraPhaseProcess, DefaultConstructor){
  std::vector<micm::Species> reactants;
  std::vector<micm::Species> products;
  micm::IntraPhaseProcess<micm::ArrheniusRateConstant> arrhenius(reactants, products, micm::ArrheniusRateConstant());
  micm::IntraPhaseProcess<micm::TroeRateConstant> troe(reactants, products, micm::TroeRateConstant());
  micm::IntraPhaseProcess<micm::ExternalRateConstant> external(reactants, products, micm::ExternalRateConstant());
}