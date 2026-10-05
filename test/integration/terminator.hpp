#include "../precision_matchers.hpp"

#include <micm/CPU.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

/// @brief A test of the "Terminator" mechanism:
///
/// Cl2 --hv--> 2 Cl
/// Cl + Cl --> Cl2
///
/// More details including analytical solution can be found here:
/// https://github.com/ESCOMP/CAM/blob/8cd44c50fe107c0b93ccd48b61eaa3d10a5b4e2f/src/chemistry/pp_terminator/chemistry.F90#L1-L434
template<class BuilderPolicy>
void TestTerminator(BuilderPolicy& builder, micm::Index number_of_grid_cells)
{
  auto cl2 = micm::Species("Cl2");
  auto cl = micm::Species("Cl");
  cl.SetProperty("absolute tolerance", micm::Real{ 1.0e-20 });
  cl2.SetProperty("absolute tolerance", micm::Real{ 1.0e-20 });

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ cl2, cl } };

  micm::Process toy_r1 = micm::ChemicalReactionBuilder()
                             .SetReactants({ cl2 })
                             .SetProducts({ micm::StoichSpecies(cl, 2.0) })
                             .SetPhase(gas_phase)
                             .SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = "toy_k1" })
                             .Build();

  constexpr micm::Real k2 = 1.0;
  micm::Process toy_r2 = micm::ChemicalReactionBuilder()
                             .SetReactants({ cl, cl })
                             .SetProducts({ micm::StoichSpecies(cl2, 1.0) })
                             .SetPhase(gas_phase)
                             .SetRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = k2 })
                             .Build();

  auto solver =
      builder.SetSystem(micm::System(gas_phase)).SetReactions(std::vector<micm::Process>{ toy_r1, toy_r2 }).Build();
  auto state = solver.GetState(number_of_grid_cells);
  // 1e-8 is below float epsilon (1.2e-7), so in a single-precision build the adaptive error norm can
  // never fall below 1 and the step size collapses to StepSizeTooSmall. Track the working precision,
  // matching the default in micm::StateParameters.
  state.SetRelativeTolerance(micm::Real{ std::is_same_v<micm::Real, double> ? 1.0e-8 : 1.0e-5 });

  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());
  std::unordered_map<std::string, std::vector<micm::Real>> concentrations{ { "Cl2", {} }, { "Cl", {} } };
  for (micm::Index i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    concentrations["Cl2"].push_back(get_double() * 1.0e-6);
    concentrations["Cl"].push_back(get_double() * 1.0e-10);
  }
  state.SetConcentrations(concentrations);

  std::unordered_map<std::string, std::vector<micm::Real>> custom_rate_constants{
    { "toy_k1", std::vector<micm::Real>(number_of_grid_cells) }
  };

  auto orig_state_vars = state.variables_;
  micm::Index steps = std::floor(2.0 * M_PI / 0.3);
  for (micm::Index step = 0; step < steps; ++step)
  {
    micm::Real lon = step * 0.3;
    state.variables_ = orig_state_vars;
    for (micm::Index i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
    {
      constexpr micm::Real k1_lat_center = M_PI * 20.0 / 180.0;
      constexpr micm::Real k1_lon_center = M_PI * 300.0 / 180.0;
      micm::Real lat = M_PI / 180.0 * (i_cell * (90.0 / number_of_grid_cells));

      micm::Real k1 = std::max<micm::Real>(
          0.0,
          std::sin(lat) * std::sin(k1_lat_center) + std::cos(lat) * std::cos(k1_lat_center) * std::cos(lon - k1_lon_center));
      custom_rate_constants["toy_k1"][i_cell] = k1;
      state.conditions_[i_cell].temperature_ = 298.0;  // K
      state.conditions_[i_cell].pressure_ = 101300.0;  // Pa
      state.conditions_[i_cell].air_density_ = 42.0;   // mol m-3
    }
    state.SetCustomRateParameters(custom_rate_constants);

    state.variables_.CopyToDevice();
    state.conditions_.CopyToDevice();
    state.custom_rate_parameters_.CopyToDevice();

    solver.UpdateStateParameters(state);

    micm::Real dt = 30.0;
    auto result = solver.Solve(dt, state);

    state.variables_.CopyToHost();

    EXPECT_EQ(result.state_, micm::SolverState::Converged);

    for (micm::Index i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
    {
      // Reference evaluated in double whatever the solver's precision: cl_f differences the nearly
      // equal quantities cl_i, det and r, and that cancellation is not resolvable in float.
      const double r = (double)custom_rate_constants["toy_k1"][i_cell] / (4.0 * (double)k2);
      const double cl_i = concentrations["Cl"][i_cell];
      const double cl2_i = concentrations["Cl2"][i_cell];
      const double cly = cl_i + 2.0 * cl2_i;
      const double det = std::sqrt(r * r + 2.0 * r * cly);
      const double e = std::exp(-4.0 * (double)k2 * det * (double)dt);
      const double l = (det * (double)k2 * (double)dt) > 1.0e-16 ? (1.0 - e) / det / (double)dt : 4.0 * (double)k2;
      const double cl_f = -l * (cl_i - det + r) * (cl_i + det + r) / (1.0 + e + (double)dt * l * (cl_i + r));
      const double cl2_f = -cl_f / 2.0;
      const double cl_ref = cl_i + (double)dt * cl_f;
      const double cl2_ref = cl2_i + (double)dt * cl2_f;
      // The 1e-6 relative term is the bound this test has always held the double build to. In float
      // it is tighter than the solver's own relative tolerance set above, so the comparison is
      // floored at the accuracy an integrated result can actually reach.
      EXPECT_REAL_SOLVE_CLOSE(
          state.variables_[i_cell][state.variable_map_["Cl"]], cl_ref, std::abs(cl_ref) * 1.0e-6 + 1.0e-14);
      EXPECT_REAL_SOLVE_CLOSE(
          state.variables_[i_cell][state.variable_map_["Cl2"]], cl2_ref, std::abs(cl2_ref) * 1.0e-6 + 1.0e-14);
    }
  }
}
