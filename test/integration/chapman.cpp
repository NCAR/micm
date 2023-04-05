#include <gtest/gtest.h>
#include <vector>

#include <micm/system/system.hpp>
#include <micm/system/phase.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/process.hpp>

namespace micm {

  struct Reaction {
    std::vector<Species> reactants_;
    std::vector<Species> products_;
    micm::RateConstant rate_constant_;
    Phase* phase_;
  };
}

TEST(SystemBuilder, DefaultConstructor){
  auto o = micm::Species("O");
  auto o1d = micm::Species("O1D");
  auto o2 = micm::Species("O2");
  auto o3 = micm::Species("O3");
  auto m = micm::Species("M");
  auto ar = micm::Species("Ar");
  auto n = micm::Species("N");
  auto h2o = micm::Species("H2O");
  auto co2 = micm::Species("CO2");

  micm::Phase gas_phase{
    std::vector<micm::Species> {
      o, o1d, o2, o3, m, ar, n, h2o, co2
    }
  };

  micm::Reaction r {
    .reactants_ = {o, o1d},
    .products_ = {h2o, m},
    .rate_constant_ = micm::ArrheniusRateConstant(),
    .phase_ = &gas_phase
  }
}