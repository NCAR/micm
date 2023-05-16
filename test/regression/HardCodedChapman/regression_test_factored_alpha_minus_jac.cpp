#include <ISO_Fortran_binding.h>

#include <gtest/gtest.h>
#include <micm/solver/chapman_ode_solver.hpp>

extern "C" {
  void Finit(void);
  void factored_alpha_minus_jac(CFI_cdesc_t * dforce_dy, double alpha, CFI_cdesc_t * LU);
}

std::vector<double> call_fortran_factored_alpha_minus_jac(std::vector<double>& dforce_dy, const double& alpha){
  std::vector<double> result{};

  CFI_CDESC_T(1) LU;
  CFI_establish((CFI_cdesc_t *)&LU, NULL,
                      CFI_attribute_pointer,
                      CFI_type_double, 0, (CFI_rank_t)1, NULL);

  CFI_CDESC_T(1) f_dforce_dy;
  CFI_index_t extent[1] = { (long int) dforce_dy.size() };
  CFI_establish((CFI_cdesc_t *)&f_dforce_dy, dforce_dy.data(),
                      CFI_attribute_other,
                      CFI_type_double, dforce_dy.size() * sizeof(double), (CFI_rank_t)1, extent);

  factored_alpha_minus_jac(
    (CFI_cdesc_t *)&f_dforce_dy, 
    alpha, 
    (CFI_cdesc_t *)&LU
  );

  for(size_t i{}; i < LU.dim[0].extent; ++i) {
    double* d = (double *) ((char *)LU.base_addr + i * LU.elem_len);
    result.push_back(*d);
  }

  return result;
}

TEST(RegressionChapmanODESolver, factored_alpha_minus_jac){
  micm::ChapmanODESolver solver{};
  std::vector<double> dforce_dy(23, 1);
  double alpha{2};

  auto LU = solver.factored_alpha_minus_jac(dforce_dy, alpha);
  auto f_LU = call_fortran_factored_alpha_minus_jac(dforce_dy, alpha);

  EXPECT_EQ(LU.size(), f_LU.size());
  for(size_t i{}; i < LU.size(); ++i) {
    EXPECT_EQ(LU[i], f_LU[i]);
  }
}