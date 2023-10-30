#include <micm/system/system.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/species.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

// based on CAMP/test/unit_rxn_data/test_rxn_surface.F90

int main(const int argc, const char* argv[])
{

  const double rxn_gamma          = 2.0e-2;  // [unitless]
  const double bar_yield          = 1.0;     // [unitless]
  const double baz_yield          = 0.4;     // [unitless]
  const double DENSITY_stuff      = 1000.0;  // [kg m-3]
  const double DENSITY_more_stuff = 1000.0;  // [kg m-3]
  const double MW_stuff           = 0.5;     // [kg mol-1]
  const double MW_more_stuff      = 0.2;     // [kg mol-1]
  const double MW_foo             = 0.04607; // [kg mol-1]
  const double Dg_foo             = 0.95e-5; // diffusion coeff [m2 s-1]
  const double sp_number          = 1.3e6;   // single-particle number concentration [# m-3]
  const double mode_GMD           = 1.0e-6;  // mode geometric mean diameter [m]
  const double mode_GSD           = 0.1;     // mode geometric standard deviation [unitless]

  return 0;
}
