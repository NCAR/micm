#include <gtest/gtest.h>

#include <micm/process/surface_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/species.hpp>
#include <micm/system/system.hpp>
#include <micm/util/constants.hpp>

template<class OdeSolverPolicy>
void test_analytical_surface_rxn(
    const std::function<OdeSolverPolicy(const micm::System&, const std::vector<micm::Process>&)> create_solver)
{
  // parameters, from CAMP/test/unit_rxn_data/test_rxn_surface.F90
  const double mode_GMD = 1.0e-6;            // mode geometric mean diameter [m]
  const double mode_GSD = 0.1;               // mode geometric standard deviation [unitless]
  const double DENSITY_stuff = 1000.0;       // [kg m-3]
  const double DENSITY_more_stuff = 1000.0;  // [kg m-3]
  const double MW_stuff = 0.5;               // [kg mol-1]
  const double MW_more_stuff = 0.2;          // [kg mol-1]
  const double MW_foo = 0.04607;             // [kg mol-1]
  const double Dg_foo = 0.95e-5;             // diffusion coefficient [m2 s-1]
  const double rxn_gamma = 2.0e-2;           // [unitless]
  const double bar_yield = 1.0;              // [unitless]
  const double baz_yield = 0.4;              // [unitless]

  // environment
  const double temperature = 272.5;  // temperature (K)
  const double pressure = 101253.3;  // pressure (Pa)

  // initial conditions
  const double conc_foo = 1.0;
  const double conc_stuff = 2.0e-3;
  const double conc_more_stuff = 3.0e-3;

  // effective radius
  double radius = mode_GMD / 2.0 * exp(5.0 * log(mode_GSD) * log(mode_GSD) / 2.0);

  // particle number concentration [# m-3]
  double number_conc = 6.0 /
                       (M_PI * std::pow(mode_GMD, 3.0) * std::exp(9.0 / 2.0 * std::log(mode_GSD) * std::log(mode_GSD))) *
                       (conc_stuff / DENSITY_stuff + conc_more_stuff / DENSITY_more_stuff);

  micm::Species foo("foo", { { "molecular weight [kg mol-1]", MW_foo }, { "diffusion coefficient [m2 s-1]", Dg_foo } });
  micm::Species bar("bar");
  micm::Species baz("baz");

  // Phase
  micm::Phase gas_phase{ std::vector<micm::Species>{ foo, bar, baz } };

  // System
  micm::System chemical_system = micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase });

  // Rate
  micm::SurfaceRateConstant surface{ { .label_ = "foo", .species_ = foo, .reaction_probability_ = rxn_gamma } };

  // Process
  micm::Process surface_process = micm::Process::create()
                                      .reactants({ foo })
                                      .products({ micm::yields(bar, bar_yield), micm::yields(baz, baz_yield) })
                                      .rate_constant(surface)
                                      .phase(gas_phase);

  auto reactions = std::vector<micm::Process>{ surface_process };

  // Solver
  // micm::RosenbrockSolver<> solver{
  //   chemical_system, reactions, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  // };

  // Solver
  OdeSolverPolicy solver = create_solver(chemical_system, reactions);

  // State
  micm::State state = solver.GetState();
  state.conditions_[0].temperature_ = temperature;
  state.conditions_[0].pressure_ = pressure;
  state.SetCustomRateParameter("foo.effective radius [m]", radius);
  state.SetCustomRateParameter("foo.particle number concentration [# m-3]", number_conc);
  state.SetConcentration(foo, conc_foo);

  // Surface reaction rate calculation
  double mean_free_speed = std::sqrt(8.0 * GAS_CONSTANT / (M_PI * MW_foo) * temperature);
  double k1 = 4.0 * number_conc * M_PI * radius * radius / (radius / Dg_foo + 4.0 / (mean_free_speed * rxn_gamma));

  double time_step = 0.1 / k1;  // s
  int nstep = 10;

  std::vector<std::vector<double>> model_conc(nstep + 1, std::vector<double>(3));
  std::vector<std::vector<double>> analytic_conc(nstep + 1, std::vector<double>(3));

  model_conc[0] = { conc_foo, 0, 0 };
  analytic_conc[0] = { conc_foo, 0, 0 };

  size_t idx_foo = 0, idx_bar = 1, idx_baz = 2;

  std::cout << std::setw(3) << "i" << std::setw(7) << "time" << std::setw(11) << "anal foo" << std::setw(11) << "model foo"
            << std::setw(11) << "anal bar" << std::setw(11) << "model bar" << std::setw(11) << "anal baz" << std::setw(11)
            << "model baz" << std::endl;

  for (int i = 1; i <= nstep; ++i)
  {
    double elapsed_solve_time = 0;

    // first iteration
    auto result = solver.Solve(time_step - elapsed_solve_time, state);
    elapsed_solve_time = result.final_time_;
    state.variables_ = result.result_;

    // further iterations
    /*
    while (elapsed_solve_time < time_step)
    {
      result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
      state.variables_ = result.result_;
    }
    */
    EXPECT_EQ(result.state_, (micm::SolverState::Converged));

    // Check surface reaction rate calculation
    EXPECT_NEAR(k1, state.rate_constants_.AsVector()[0], 1e-8);

    model_conc[i] = result.result_.AsVector();

    double time = i * time_step;
    analytic_conc[i][idx_foo] = conc_foo * std::exp(-k1 * time);
    analytic_conc[i][idx_bar] = bar_yield * (1.0 - analytic_conc[i][idx_foo]);
    analytic_conc[i][idx_baz] = baz_yield * (1.0 - analytic_conc[i][idx_foo]);

    // Check concentrations
    EXPECT_NEAR(analytic_conc[i][idx_foo], model_conc[i][idx_foo], 1e-5);
    EXPECT_NEAR(analytic_conc[i][idx_bar], model_conc[i][idx_bar], 1e-5);
    EXPECT_NEAR(analytic_conc[i][idx_baz], model_conc[i][idx_baz], 1e-5);

    std::cout << std::setw(3) << i << "  " << std::fixed << std::setprecision(2) << std::setw(5) << time << "  "
              << std::fixed << std::setprecision(7) << analytic_conc[i][idx_foo] << "  " << model_conc[i][idx_foo] << "  "
              << analytic_conc[i][idx_bar] << "  " << model_conc[i][idx_bar] << "  " << analytic_conc[i][idx_baz] << "  "
              << model_conc[i][idx_baz] << "  " << std::endl;
  }
}
