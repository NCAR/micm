#include <ISO_Fortran_binding.h>

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/chapman_ode_solver.hpp>

#include <gtest/gtest.h>

extern "C" {
  double arrhenius_rate(
    double temperature, 
    double pressure, 
    double a, 
    double b, 
    double c, 
    double d, 
    double e
  );

  void reaction_rate_constants(
    double temperature, 
    double pressure, 
    double photorate_1, 
    double photorate_2, 
    double photorate_3, 
    CFI_cdesc_t * new_number_densities
  );
}

double call_fortran_rate(double temperature, double pressure, micm::ArrheniusRateConstantParameters params){
  return arrhenius_rate(
    temperature, pressure,
    params.A_, params.B_, params.C_, params.D_, params.E_
  );
}

std::vector<double> call_reaction_rate_constants(double temperature, double pressure, double photorate_1, double photorate_2, double photorate_3){
  std::vector<double> result{};

  CFI_CDESC_T(1) f_reaction_rate_constants;
  CFI_establish((CFI_cdesc_t *)&f_reaction_rate_constants, NULL,
                      CFI_attribute_pointer,
                      CFI_type_double, 0, (CFI_rank_t)1, NULL);

  reaction_rate_constants(
    temperature, pressure, photorate_1, photorate_2, photorate_3,
    (CFI_cdesc_t *)&f_reaction_rate_constants
  );

  for(size_t i{}; i < f_reaction_rate_constants.dim[0].extent; ++i) {
    double* d = (double *) ((char *)f_reaction_rate_constants.base_addr + i * f_reaction_rate_constants.elem_len);
    result.push_back(*d);
  }

  return result;
}

TEST(ArrheniusRateRegressionTest, AllZeros){
  micm::State state {0,0};
  state.temperature_ = 301.24; // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::ArrheniusRateConstant zero{};
  auto k = zero.calculate(state, params);
  auto k_fortran = call_fortran_rate(1, 1, micm::ArrheniusRateConstantParameters());
  EXPECT_NEAR(k, k_fortran, 1e-10);

}

TEST(ArrheniusRateRegressionTest, A1RestZero){
  micm::State state {0,0};
  state.temperature_ = 301.24; // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 1;
  micm::ArrheniusRateConstant basic(parameters);
  auto k = basic.calculate(state, params);
  auto k_fortran = call_fortran_rate(1, 1, parameters);
  EXPECT_NEAR(k, k_fortran, 1e-10);
}

// values from https://jpldataeval.jpl.nasa.gov/pdf/JPL_00-03.pdf
TEST(ArrheniusRateRegressionTest, O1D){
  micm::State state {0,0};
  state.temperature_ = 301.24; // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 2.2e-10;
  micm::ArrheniusRateConstant o1d(parameters);
  auto k = o1d.calculate(state, params);
  auto k_fortran = call_fortran_rate(298, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 1e-10);
}

// O + HO2 -> OH + O2
TEST(ArrheniusRateRegressionTest, HOx){
  micm::State state {0,0};
  state.temperature_ = 301.24; // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 3e-11;
  parameters.C_ = -200;
  micm::ArrheniusRateConstant hox(parameters);
  auto k = hox.calculate(state, params);
  auto k_fortran = call_fortran_rate(301.24, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 1e-10);
}

// OH + HCl → H2O + Cl
TEST(ArrheniusRateRegressionTest, ClOx){
  micm::State state {0,0};
  state.temperature_ = 301.24; // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 2.6e-12;
  parameters.C_ = -350;
  micm::ArrheniusRateConstant clox(parameters);
  auto k = clox.calculate(state, params);
  auto k_fortran = call_fortran_rate(301.24, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 1e-10);
}

// k_N2_O1D_1: N2 + O1D -> 1*O + 1*N2
TEST(ArrheniusRateRegressionTest, N2_O1D_1){
  micm::State state {0,0};
  state.temperature_ = 301.24; // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 2.15e-11;
  parameters.C_ = 110;
  micm::ArrheniusRateConstant rate(parameters);
  auto k = rate.calculate(state, params);
  auto k_fortran = call_fortran_rate(301.24, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 1e-10);
}

// k_O1D_O2_1: O1D + O2 -> 1*O + 1*O2
TEST(ArrheniusRateRegressionTest, O1D_O2_1){
  micm::State state {0,0};
  state.temperature_ = 301.24; // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 3.3e-11;
  parameters.C_ = 55;
  micm::ArrheniusRateConstant rate(parameters);
  auto k = rate.calculate(state, params);
  auto k_fortran = call_fortran_rate(301.24, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 1e-10);
}

// k_O_O3_1: O + O3 -> 2*O2
TEST(ArrheniusRateRegressionTest, O_O3_1){
  micm::State state {0,0};
  state.temperature_ = 301.24; // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 8e-12;
  parameters.C_ = -2060;
  micm::ArrheniusRateConstant rate(parameters);
  auto k = rate.calculate(state, params);
  auto k_fortran = call_fortran_rate(301.24, 100'000, parameters);
  EXPECT_NEAR(k, k_fortran, 1e-10);
}

// k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
TEST(ArrheniusRateRegressionTest, M_O_O2_1){
  micm::State state {0,0};
  state.temperature_ = 273.14; // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::ArrheniusRateConstantParameters parameters;
  parameters.A_ = 6e-34;
  parameters.B_ = -2.4;
  micm::ArrheniusRateConstant rate(parameters);
  auto k = rate.calculate(state, params);
  auto k_fortran = call_fortran_rate(273.14, 100'000, parameters);
  EXPECT_EQ(k, k_fortran);
}

TEST(ArrheniusRateRegressionTest, integration){
  micm::State state {0,0};
  state.temperature_ = 273.15; // [K]
  std::vector<double>::const_iterator params = state.custom_rate_parameters_.begin();
  micm::ChapmanODESolver solver{};
  double temperature = 273.15;
  double pressure = 1000 * 100; // 1000 hPa

  solver.calculate_rate_constants(state);
  auto f_reaction_rate_constants = call_reaction_rate_constants(temperature, pressure, solver.rate_constants_[0], solver.rate_constants_[1], solver.rate_constants_[2]);

  EXPECT_NEAR(solver.rate_constants_[0], f_reaction_rate_constants[0], 1e-8);
  EXPECT_NEAR(solver.rate_constants_[1], f_reaction_rate_constants[1], 1e-8);
  EXPECT_NEAR(solver.rate_constants_[2], f_reaction_rate_constants[2], 1e-8);
  EXPECT_NEAR(solver.rate_constants_[3], f_reaction_rate_constants[3], 1e-8);
  EXPECT_NEAR(solver.rate_constants_[4], f_reaction_rate_constants[4], 1e-8);
  EXPECT_NEAR(solver.rate_constants_[5], f_reaction_rate_constants[5], 1e-8);
  EXPECT_NEAR(solver.rate_constants_[6], f_reaction_rate_constants[6], 1e-8);
}
