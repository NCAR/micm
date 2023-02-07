#include <ISO_Fortran_binding.h>

#include <gtest/gtest.h>
#include <micm/solver/chapman_ode_solver.hpp>

extern "C" {
  void solve(CFI_cdesc_t * jacobian, CFI_cdesc_t * b, CFI_cdesc_t * x);
}

std::vector<double> call_fortran_solve(std::vector<double>& jacobian, std::vector<double>& b){
  std::vector<double> result{};

  CFI_CDESC_T(1) f_x;
  CFI_establish((CFI_cdesc_t *)&f_x, NULL,
                      CFI_attribute_pointer,
                      CFI_type_double, 0, (CFI_rank_t)1, NULL);

  CFI_CDESC_T(1) f_jacobian;
  CFI_index_t extent[1] = { (long int) jacobian.size() };
  CFI_establish((CFI_cdesc_t *)&f_jacobian, jacobian.data(),
                      CFI_attribute_other,
                      CFI_type_double, jacobian.size() * sizeof(double), (CFI_rank_t)1, extent);

  CFI_CDESC_T(1) f_b;
  CFI_index_t b_extent[1] = { (long int) b.size() };
  CFI_establish((CFI_cdesc_t *)&f_b, b.data(),
                      CFI_attribute_other,
                      CFI_type_double, b.size() * sizeof(double), (CFI_rank_t)1, b_extent);

  solve(
    (CFI_cdesc_t *)&f_jacobian, 
    (CFI_cdesc_t *)&f_b,
    (CFI_cdesc_t *)&f_x
  );

  for(size_t i{}; i < f_x.dim[0].extent; ++i) {
    double* d = (double *) ((char *)f_x.base_addr + i * f_x.elem_len);
    result.push_back(*d);
  }

  return result;
}

TEST(RegressionChapmanODESolver, lin_solve){
  micm::ChapmanODESolver solver{};
  std::vector<double> jacobian(23, 1), b(23, 0.5);

  auto solved = solver.lin_solve(b, jacobian);
  auto f_solved = call_fortran_solve(jacobian, b);

  EXPECT_EQ(solved.size(), f_solved.size());
  for(size_t i{}; i < solved.size(); ++i) {
    EXPECT_EQ(solved[i], f_solved[i]);
  }
}