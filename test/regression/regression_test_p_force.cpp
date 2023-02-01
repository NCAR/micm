#include <ISO_Fortran_binding.h>

#include <gtest/gtest.h>
#include <micm/solver/chapman_ode_solver.hpp>

extern "C" {
  void Finit(void);
  void p_force( CFI_cdesc_t * );
}

TEST(RegressionChapmanODESolver, PForce){
  micm::ChapmanODESolver solver{};
  std::vector<double> rate_constants(9, 0);
  std::vector<double> number_densities(9, 0);
  double number_density_air{};

  auto forcing = solver.p_force(rate_constants, number_densities, number_density_air);

  EXPECT_EQ(forcing[0], 0);
  EXPECT_EQ(forcing[1], 0);
  EXPECT_EQ(forcing[2], 0);
  EXPECT_EQ(forcing[3], 0);
  EXPECT_EQ(forcing[4], 0);
}