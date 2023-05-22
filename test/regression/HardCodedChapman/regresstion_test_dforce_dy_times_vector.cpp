#include <ISO_Fortran_binding.h>
#include <gtest/gtest.h>

#include <micm/solver/chapman_ode_solver.hpp>

extern "C"
{
  void dforce_dy_times_vector(CFI_cdesc_t *dforce_dy, CFI_cdesc_t *vector, CFI_cdesc_t *product);
}

std::vector<double> call_fortran_dforce_dy_times_vector(std::vector<double> &dforce_dy, std::vector<double> &vector)
{
  std::vector<double> result{};

  CFI_CDESC_T(1) product;
  CFI_establish((CFI_cdesc_t *)&product, NULL, CFI_attribute_pointer, CFI_type_double, 0, (CFI_rank_t)1, NULL);

  CFI_CDESC_T(1) f_dforce_dy;
  CFI_index_t extent[1] = { (long int)dforce_dy.size() };
  CFI_establish(
      (CFI_cdesc_t *)&f_dforce_dy,
      dforce_dy.data(),
      CFI_attribute_other,
      CFI_type_double,
      dforce_dy.size() * sizeof(double),
      (CFI_rank_t)1,
      extent);

  CFI_CDESC_T(1) f_vector;
  CFI_index_t vector_extent[1] = { (long int)vector.size() };
  CFI_establish(
      (CFI_cdesc_t *)&f_vector,
      vector.data(),
      CFI_attribute_other,
      CFI_type_double,
      vector.size() * sizeof(double),
      (CFI_rank_t)1,
      vector_extent);

  dforce_dy_times_vector((CFI_cdesc_t *)&f_dforce_dy, (CFI_cdesc_t *)&f_vector, (CFI_cdesc_t *)&product);

  for (size_t i{}; i < product.dim[0].extent; ++i)
  {
    double *d = (double *)((char *)product.base_addr + i * product.elem_len);
    result.push_back(*d);
  }

  return result;
}

TEST(RegressionChapmanODESolver, factored_alpha_minus_jac)
{
  micm::ChapmanODESolver solver{};
  std::vector<double> dforce_dy(23, 1);
  std::vector<double> vector(23, 0.5);

  auto product = solver.dforce_dy_times_vector(dforce_dy, vector);
  auto f_product = call_fortran_dforce_dy_times_vector(dforce_dy, vector);

  EXPECT_EQ(product.size(), f_product.size());
  for (size_t i{}; i < product.size(); ++i)
  {
    EXPECT_EQ(product[i], f_product[i]);
  }
}