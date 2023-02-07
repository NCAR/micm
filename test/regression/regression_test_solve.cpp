#include <ISO_Fortran_binding.h>

#include <gtest/gtest.h>
#include <micm/solver/chapman_ode_solver.hpp>

extern "C" {
  void solve();
}

std::vector<double> call_fortran_solve(){
  solve();

  return std::vector<double>();
}

TEST(RegressionChapmanODESolver, lin_solve){
  micm::ChapmanODESolver solver{};
  std::vector<double> number_densities(23, 5e-8);
  double number_density_air = 2.7e19;

  solver.calculate_rate_constants(273.15, 100000);

  auto results = solver.Solve(0, 0.1, number_densities, number_density_air);
  auto f_results = call_fortran_solve();

  std::cout << "solver state: " << micm::state_to_string(results.state_) << "\n";
}