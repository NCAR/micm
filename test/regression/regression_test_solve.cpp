#include <ISO_Fortran_binding.h>

#include <gtest/gtest.h>
#include <micm/solver/chapman_ode_solver.hpp>

extern "C" {
  void solve(
    double temperature, double pressure, double number_density_air, 
    double time_start, double time_end, 
    CFI_cdesc_t * number_densities,
    CFI_cdesc_t * new_number_densities
  );
}

std::vector<double> call_fortran_solve(
    double temperature, double pressure, double number_density_air, 
    double time_start, double time_end, 
    std::vector<double> number_densities)
  {

  std::vector<double> result{};

  CFI_CDESC_T(1) f_number_densities;
  CFI_index_t extent[1] = { (long int) number_densities.size() };
  CFI_establish((CFI_cdesc_t *)&f_number_densities, number_densities.data(),
                      CFI_attribute_other,
                      CFI_type_double, number_densities.size() * sizeof(double), (CFI_rank_t)1, extent);

  CFI_CDESC_T(1) f_new_number_densities;
  CFI_establish((CFI_cdesc_t *)&f_new_number_densities, NULL,
                      CFI_attribute_pointer,
                      CFI_type_double, 0, (CFI_rank_t)1, NULL);
  solve(
    temperature, pressure, number_density_air, time_start, time_end, 
    (CFI_cdesc_t *)&f_number_densities,
    (CFI_cdesc_t *)&f_new_number_densities
  );

  for(size_t i{}; i < f_new_number_densities.dim[0].extent; ++i) {
    double* d = (double *) ((char *)f_new_number_densities.base_addr + i * f_new_number_densities.elem_len);
    result.push_back(*d);
  }

  return result;
}

TEST(RegressionChapmanODESolver, solve){
  micm::ChapmanODESolver solver{};
  std::vector<double> number_densities = { 1,    3.92e-1, 1.69e-2, 0,     3.29e1, 0,     0,   8.84, 0};
                                         //"M"   "Ar"     "CO2",   "H2O", "N2",   "O1D", "O", "O2", "O3",
  double number_density_air = 2.7e19;
  double temperature = 273.15;
  double pressure = 1000 * 100; // 1000 hPa
  double time_start = 0;
  double time_end = 1;

  solver.calculate_rate_constants(temperature, pressure);

  auto f_results = call_fortran_solve(temperature, pressure, number_density_air, time_start, time_end, number_densities);
  auto results = solver.Solve(time_start, time_end, number_densities, number_density_air);

  for(const auto& elem : f_results){
    std::cout << elem << " ";
  }
  std::cout << "\n";

  // EXPECT_DEATH(
  //   call_fortran_solve(temperature, pressure, number_density_air, time_start, time_end, number_densities),
  //   "asdf"
  // );

  std::cout << "solver state: " << micm::state_to_string(results.state_) << "\n";
}