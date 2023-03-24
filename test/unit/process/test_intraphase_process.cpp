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

TEST(IntraPhaseProcess, ConstructorWithArguments){
  std::vector<micm::Species> reactants;
  std::vector<micm::Species> products;
  micm::IntraPhaseProcess<micm::ArrheniusRateConstant> arrhenius(reactants, products, micm::ArrheniusRateConstant());
  micm::IntraPhaseProcess<micm::TroeRateConstant> troe(reactants, products, micm::TroeRateConstant());
  micm::IntraPhaseProcess<micm::ExternalRateConstant> external(reactants, products, micm::ExternalRateConstant());
}

TEST(IntraPhaseProcess, CopyConstructor) {
  std::vector<micm::Species> reactants;
  std::vector<micm::Species> products;
  micm::IntraPhaseProcess<micm::ArrheniusRateConstant> arrhenius(reactants, products, micm::ArrheniusRateConstant());
  micm::IntraPhaseProcess<micm::TroeRateConstant> troe(reactants, products, micm::TroeRateConstant());
  micm::IntraPhaseProcess<micm::ExternalRateConstant> external(reactants, products, micm::ExternalRateConstant());

  auto arrhenius2(arrhenius);
  auto troe2(troe);
  auto external2(external);
}

TEST(IntraPhaseProcess, MoveConstructor) {
  std::vector<micm::Species> reactants;
  std::vector<micm::Species> products;

  auto arrhenius(micm::IntraPhaseProcess<micm::ArrheniusRateConstant>(reactants, products, micm::ArrheniusRateConstant()));
  auto troe(micm::IntraPhaseProcess<micm::TroeRateConstant>(reactants, products, micm::TroeRateConstant()));
  auto external(micm::IntraPhaseProcess<micm::ExternalRateConstant>(reactants, products, micm::ExternalRateConstant()));
}

TEST(IntraPhaseProcess, CopyAssignment) {
  std::vector<micm::Species> reactants;
  std::vector<micm::Species> products;
  micm::IntraPhaseProcess<micm::ArrheniusRateConstant> arrhenius(reactants, products, micm::ArrheniusRateConstant());
  micm::IntraPhaseProcess<micm::TroeRateConstant> troe(reactants, products, micm::TroeRateConstant());
  micm::IntraPhaseProcess<micm::ExternalRateConstant> external(reactants, products, micm::ExternalRateConstant());

  auto arrhenius2 = arrhenius;
  auto troe2 = troe;
  auto external2 = external;
}

TEST(IntraPhaseProcess, MoveAssignment) {
  std::vector<micm::Species> reactants;
  std::vector<micm::Species> products;

  auto arrhenius = micm::IntraPhaseProcess<micm::ArrheniusRateConstant>(reactants, products, micm::ArrheniusRateConstant());
  auto troe = micm::IntraPhaseProcess<micm::TroeRateConstant>(reactants, products, micm::TroeRateConstant());
  auto external = micm::IntraPhaseProcess<micm::ExternalRateConstant>(reactants, products, micm::ExternalRateConstant());
}