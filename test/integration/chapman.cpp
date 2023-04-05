#include <gtest/gtest.h>
#include <vector>
#include <pair>

#include <micm/system/system.hpp>
#include <micm/system/phase.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/process.hpp>

TEST(SystemBuilder, DefaultConstructor){
  auto o = micm::Species("O");
  auto o1d = micm::Species("O1D");
  auto o2 = micm::Species("O2");
  auto o3 = micm::Species("O3");
  auto m = micm::Species("M");
  auto ar = micm::Species("Ar");
  auto n2 = micm::Species("N2");
  auto h2o = micm::Species("H2O");
  auto co2 = micm::Species("CO2");

  micm::Phase gas_phase{
    std::vector<micm::Species> {
      o, o1d, o2, o3, m, ar, n, h2o, co2
    }
  };

  micm::Process r1 {
    .reactants_ = {o1d, n2},
    .products_ = { yields(o, 1), yields(n2, 1) },
    .rate_constant_ = micm::ArrheniusRateConstant(
      micm::ArrheniusRateConstantParameters { .A_ = 2.15e-11, .C_=110 }
    ),
    .phase_ = &gas_phase
  };

  micm::Process r2 {
    .reactants_ = {o1d, o2},
    .products_ = { yields(o, 1), yields(o2, 1) },
    .rate_constant_ = micm::ArrheniusRateConstant(
      micm::ArrheniusRateConstantParameters { .A_ = 3.3e-11, .C_=55 }
    ),
    .phase_ = &gas_phase
  };

  micm::Process r3 {
    .reactants_ = {o, o3},
    .products_ = { yields(o2, 2) },
    .rate_constant_ = micm::ArrheniusRateConstant(
      micm::ArrheniusRateConstantParameters { .A_ = 8e-12, .C_=-2060 }
    ),
    .phase_ = &gas_phase
  };

  micm::Process r4 {
    .reactants_ = {o, o2, M},
    .products_ = { yields(o3, 1), yields(m, 1) },
    .rate_constant_ = micm::ArrheniusRateConstant(
      micm::ArrheniusRateConstantParameters { .A_ = 6.0e-34, .C_=2.4 }
    ),
    .phase_ = &gas_phase
  };

  using std::pair<micm::Species, double = 1> = yields;
  micm::Process photo_1 {
    .reactants_ = { o2 },
    .products_ = { yields(o, 2) },
    .rate_constant_ = micm::PhotolysisRateConstant(),
    .phase_ = &gas_phase
  };

  micm::Process photo_2 {
    .reactants_ = { o3 },
    .products_ = { yields(o1d, 1), yields(o2, 1) },
    .rate_constant_ = micm::PhotolysisRateConstant(),
    .phase_ = &gas_phase
  };

  micm::Process photo_3 {
    .reactants_ = { o3 },
    .products_ = { yields(o, 1), yields(o2, 1) },
    .rate_constant_ = micm::PhotolysisRateConstant(),
    .phase_ = &gas_phase
  };

  struct State  {
    double temperature;
    double pressure; 
    std::vector<double> concentrations;

    std::function<> get_jacobian();
    std::function<> get_forcing();
  };

  System system{gas_phase};

  State state{system, r1, ..., photo_3 };

  state.concentrations = {1, 2 };
  state.temperature = 2;
  state.pressure = 3;
  state.photo_rates = {1, 4};

  Rosenbrock solver( state );

  for(int t {}, t < 100; ++t)
  {
    state.update_photo_rates();
    solver.solve( state );
    // output state
  }
}

class MICM {
  Solver* solver_;
  State* state_;
}

// assumes that photo_rates, matches order in c++ already
void fortran_solve(void* micm_address, double* concentrations, double temperature, double pressure, double[] photo_rates) {
  MICM* micm = std::reinterpret_pointer_cast<MICM>(micm_address);

  micm->state_->photo_rates = photo_rates;
  micm->state_->concentrations = concentrations;
  //temp and pres

  micm->solver_->solve(state);
}