#include <gtest/gtest.h>
#include <math.h>

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/system.hpp>
#include <random>
#include <utility>
#include <vector>

/// @brief A test of the "Terminator" mechanism:
///
/// Cl2 --hv--> 2 Cl
/// Cl + Cl --> Cl2
///
/// More details including analytical solution can be found here:
/// https://github.com/ESCOMP/CAM/blob/8cd44c50fe107c0b93ccd48b61eaa3d10a5b4e2f/src/chemistry/pp_terminator/chemistry.F90#L1-L434
template<template<class> class MatrixPolicy, class OdeSolverPolicy>
void TestTerminator(
    const std::function<OdeSolverPolicy(const micm::System&, const std::vector<micm::Process>&)> create_solver,
    std::size_t number_of_grid_cells)
{
  auto cl2 = micm::Species("Cl2");
  auto cl = micm::Species("Cl");

  micm::Phase gas_phase{ std::vector<micm::Species>{ cl2, cl } };

  micm::Process toy_r1 = micm::Process::create()
                             .reactants({ cl2 })
                             .products({ micm::Yield(cl, 2.0) })
                             .phase(gas_phase)
                             .rate_constant(micm::UserDefinedRateConstant({ .label_ = "toy_k1" }));

  constexpr double k2 = 1.0;
  micm::Process toy_r2 = micm::Process::create()
                             .reactants({ cl, cl })
                             .products({ micm::Yield(cl2, 1.0) })
                             .phase(gas_phase)
                             .rate_constant(micm::ArrheniusRateConstant({ .A_ = k2 }));

  auto solver = create_solver(
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }), std::vector<micm::Process>{ toy_r1, toy_r2 });
  auto state = solver.GetState();

  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());
  std::unordered_map<std::string, std::vector<double>> concentrations{ { "Cl2", {} }, { "Cl", {} } };
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    concentrations["Cl2"].push_back(get_double() * 1.0e-6);
    concentrations["Cl"].push_back(get_double() * 1.0e-10);
  }
  state.SetConcentrations(concentrations);

  std::unordered_map<std::string, std::vector<double>> custom_rate_constants{
    { "toy_k1", std::vector<double>(number_of_grid_cells) }
  };
  for (double lon = 0.0; lon < 2.0 * M_PI; lon += 0.3)
  {
    for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
    {
      constexpr double k1_lat_center = M_PI * 20.0 / 180.0;
      constexpr double k1_lon_center = M_PI * 300.0 / 180.0;
      double lat = M_PI / 180.0 * (i_cell * (90.0 / number_of_grid_cells));

      double k1 = std::max(
          0.0,
          std::sin(lat) * std::sin(k1_lat_center) + std::cos(lat) * std::cos(k1_lat_center) * std::cos(lon - k1_lon_center));
      custom_rate_constants["toy_k1"][i_cell] = k1;
      state.conditions_[i_cell].temperature_ = 298.0;  // K
      state.conditions_[i_cell].pressure_ = 101300.0;  // Pa
      state.conditions_[i_cell].air_density_ = 42.0;   // mol m-3
    }
    state.SetCustomRateParameters(custom_rate_constants);

    double dt = 30.0;
    auto result = solver.Solve(dt, state);

    EXPECT_EQ(result.state_, micm::SolverState::Converged);

    for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
    {
      double r = custom_rate_constants["toy_k1"][i_cell] / (4.0 * k2);
      double cl_i = concentrations["Cl"][i_cell];
      double cl2_i = concentrations["Cl2"][i_cell];
      double cly = cl_i + 2.0 * cl2_i;
      double det = std::sqrt(r * r + 2.0 * r * cly);
      double e = std::exp(-4.0 * k2 * det * dt);
      double l = (det * k2 * dt) > 1.0e-16 ? (1.0 - e) / det / dt : 4.0 * k2;
      double cl_f = -l * (cl_i - det + r) * (cl_i + det + r) / (1.0 + e + dt * l * (cl_i + r));
      double cl2_f = -cl_f / 2.0;
      EXPECT_NEAR(
          result.result_[i_cell][state.variable_map_["Cl"]], cl_i + dt * cl_f, (cl_i + dt * cl_f) * 1.0e-7 + 1.0e-14);
      EXPECT_NEAR(
          result.result_[i_cell][state.variable_map_["Cl2"]], cl2_i + dt * cl2_f, (cl2_i + dt * cl2_f) * 1.0e-7 + 1.0e-14);
    }
  }
}
