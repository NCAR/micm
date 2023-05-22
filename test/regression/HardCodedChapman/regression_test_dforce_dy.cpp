#include <ISO_Fortran_binding.h>
#include <gtest/gtest.h>

#include <micm/solver/chapman_ode_solver.hpp>

extern "C"
{
  void dforce_dy(
      CFI_cdesc_t *rate_constants,
      CFI_cdesc_t *number_densities,
      double number_density_air,
      CFI_cdesc_t *new_number_densities);
}

std::vector<double>
call_dforce_dy(std::vector<double> rate_constants, std::vector<double> number_densities, double number_density_air)
{
  std::vector<double> result{};

  CFI_CDESC_T(1) f_rate_constants;
  CFI_index_t extent[1] = { (long int)rate_constants.size() };
  CFI_establish(
      (CFI_cdesc_t *)&f_rate_constants,
      rate_constants.data(),
      CFI_attribute_other,
      CFI_type_double,
      rate_constants.size() * sizeof(double),
      (CFI_rank_t)1,
      extent);

  CFI_CDESC_T(1) f_number_densities;
  CFI_index_t extent_num[1] = { (long int)number_densities.size() };
  CFI_establish(
      (CFI_cdesc_t *)&f_number_densities,
      number_densities.data(),
      CFI_attribute_other,
      CFI_type_double,
      number_densities.size() * sizeof(double),
      (CFI_rank_t)1,
      extent_num);

  CFI_CDESC_T(1) f_dforce_dy;
  CFI_establish((CFI_cdesc_t *)&f_dforce_dy, NULL, CFI_attribute_pointer, CFI_type_double, 0, (CFI_rank_t)1, NULL);
  dforce_dy(
      (CFI_cdesc_t *)&f_rate_constants, (CFI_cdesc_t *)&f_number_densities, number_density_air, (CFI_cdesc_t *)&f_dforce_dy);

  for (size_t i{}; i < f_dforce_dy.dim[0].extent; ++i)
  {
    double *d = (double *)((char *)f_dforce_dy.base_addr + i * f_dforce_dy.elem_len);
    result.push_back(*d);
  }

  return result;
}

TEST(RegressionChapmanODESolver, solve)
{
  micm::ChapmanODESolver solver{};
  std::vector<double> number_densities = { 1, 3.92e-1, 1.69e-2, 0, 3.29e1, 0, 0, 8.84, 0 };
  //"M"   "Ar"     "CO2",   "H2O", "N2",   "O1D", "O", "O2", "O3",
  std::vector<double> rate_constants(number_densities.size(), 5e-7);
  double number_density_air = 2.7e19;

  auto results = solver.dforce_dy(rate_constants, number_densities, number_density_air);
  auto f_results = call_dforce_dy(rate_constants, number_densities, number_density_air);

  EXPECT_EQ(results.size(), f_results.size());
  EXPECT_EQ(results[0], f_results[0]);
  EXPECT_EQ(results[1], f_results[1]);
  EXPECT_EQ(results[2], f_results[2]);
  EXPECT_EQ(results[3], f_results[3]);
  EXPECT_EQ(results[4], f_results[4]);
  EXPECT_EQ(results[5], f_results[5]);
  EXPECT_EQ(results[6], f_results[6]);
  EXPECT_EQ(results[7], f_results[7]);
  EXPECT_EQ(results[8], f_results[8]);
  EXPECT_EQ(results[9], f_results[9]);
  EXPECT_EQ(results[10], f_results[10]);
  EXPECT_EQ(results[11], f_results[11]);
  EXPECT_EQ(results[12], f_results[12]);
  EXPECT_EQ(results[13], f_results[13]);
  EXPECT_EQ(results[14], f_results[14]);
  EXPECT_EQ(results[15], f_results[15]);
  EXPECT_EQ(results[16], f_results[16]);
  EXPECT_EQ(results[17], f_results[17]);
  EXPECT_EQ(results[18], f_results[18]);
  EXPECT_EQ(results[19], f_results[19]);
  EXPECT_EQ(results[20], f_results[20]);
  EXPECT_EQ(results[21], f_results[21]);
  EXPECT_EQ(results[22], f_results[22]);
}