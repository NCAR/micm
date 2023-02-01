#include <ISO_Fortran_binding.h>

#include <gtest/gtest.h>
#include <micm/solver/chapman_ode_solver.hpp>

extern "C" {
  void Finit(void);
  void p_force(CFI_cdesc_t * rate_constants, CFI_cdesc_t * number_densities, double number_density_air, CFI_cdesc_t * force);
}

std::vector<double> call_fortran_p_force(std::vector<double>& rate_constants, std::vector<double>& number_densities, const double& number_density_air){
  std::vector<double> result{};

  CFI_CDESC_T(1) f_p_force;
  CFI_establish((CFI_cdesc_t *)&f_p_force, NULL,
                      CFI_attribute_pointer,
                      CFI_type_double, 0, (CFI_rank_t)1, NULL);

  CFI_CDESC_T(1) f_rate_constants;
  CFI_index_t rate_constants_extent[1] = { (long int) rate_constants.size() };
  CFI_establish((CFI_cdesc_t *)&f_rate_constants, rate_constants.data(),
                      CFI_attribute_other,
                      CFI_type_double, rate_constants.size() * sizeof(double), (CFI_rank_t)1, rate_constants_extent);

  CFI_CDESC_T(1) f_number_densities;
  CFI_index_t number_densities_extent[1] = { (long int) number_densities.size() };
  CFI_establish((CFI_cdesc_t *)&f_number_densities, number_densities.data(),
                      CFI_attribute_other,
                      CFI_type_double, number_densities.size() * sizeof(double), (CFI_rank_t)1, number_densities_extent);

  p_force(
    (CFI_cdesc_t *)&f_rate_constants, 
    (CFI_cdesc_t *)&f_number_densities, 
    number_density_air, 
    (CFI_cdesc_t *)&f_p_force
  );

  for(size_t i{}; i < f_p_force.dim[0].extent; ++i) {
    double* d = (double *) ((char *)f_p_force.base_addr + i * f_p_force.elem_len);
    result.push_back(*d);
  }

  return result;
}

TEST(RegressionChapmanODESolver, simple_p_force){
  micm::ChapmanODESolver solver{};
  std::vector<double> rate_constants(9, 1);
  std::vector<double> number_densities(9, 1);
  double number_density_air{};

  auto forcing = solver.p_force(rate_constants, number_densities, number_density_air);
  auto f_forcing = call_fortran_p_force(rate_constants, number_densities, number_density_air);

  EXPECT_EQ(forcing.size(), f_forcing.size());
  for(size_t i{}; i < forcing.size(); ++i) {
    EXPECT_EQ(forcing[i], f_forcing[i]);
  }
}

TEST(RegressionChapmanODESolver, smaller_p_force){
  micm::ChapmanODESolver solver{};
  std::vector<double> rate_constants(9, 3e-8);
  std::vector<double> number_densities(9, 5e-6);
  double number_density_air{6e-14};

  auto forcing = solver.p_force(rate_constants, number_densities, number_density_air);
  auto f_forcing = call_fortran_p_force(rate_constants, number_densities, number_density_air);

  EXPECT_EQ(forcing.size(), f_forcing.size());
  for(size_t i{}; i < forcing.size(); ++i) {
    EXPECT_EQ(forcing[i], f_forcing[i]);
  }
}