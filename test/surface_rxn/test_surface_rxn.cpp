#include <micm/system/system.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/species.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

int main(const int argc, const char* argv[])
{
  // parameters, from CAMP/test/unit_rxn_data/test_rxn_surface.F90
  const double mode_GMD           = 1.0e-6;  // mode geometric mean diameter [m]
  const double mode_GSD           = 0.1;     // mode geometric standard deviation [unitless]
  const double DENSITY_stuff      = 1000.0;  // [kg m-3]
  const double DENSITY_more_stuff = 1000.0;  // [kg m-3]
  const double MW_stuff           = 0.5;     // [kg mol-1]
  const double MW_more_stuff      = 0.2;     // [kg mol-1]
  const double MW_foo             = 0.04607; // [kg mol-1]
  const double Dg_foo             = 0.95e-5; // diffusion coefficient [m2 s-1]
  const double rxn_gamma          = 2.0e-2;  // [unitless]
  const double bar_yield          = 1.0;     // [unitless]
  const double baz_yield          = 0.4;     // [unitless]

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
  double number_conc = 6.0 / (M_PI * std::pow(mode_GMD, 3.0)
    * std::exp(9.0/2.0 * std::log(mode_GSD) * std::log(mode_GSD) ))
    * (conc_stuff / DENSITY_stuff + conc_more_stuff / DENSITY_more_stuff);

  micm::Species foo("foo",
    { { "molecular weight [kg mol-1]", MW_foo },
    { "diffusion coefficient [m2 s-1]", Dg_foo } });
  micm::Species bar("bar");
  micm::Species baz("baz");

  auto state_parameters_ = micm::StateParameters{
    .number_of_grid_cells_ = 1,
    .number_of_rate_constants_ = 1,
    .variable_names_ = { "surface" },
    .custom_rate_parameter_labels_
      = { "effective radius [m]",
          "particle number concentration [# m-3]" },
  };

  micm::State state{ state_parameters_ };
  state.custom_rate_parameters_[0][0] = radius;
  state.custom_rate_parameters_[0][1] = number_conc;
  state.conditions_[0].temperature_ = temperature;
  std::vector<double>::const_iterator params = state.custom_rate_parameters_[0].begin();
  micm::SurfaceRateConstant surface{
   { .label_ = "foo", .species_ = foo, .reaction_probability_ = rxn_gamma } };

  return 0;
}
