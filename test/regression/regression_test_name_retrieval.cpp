#include <string_view>
#include <vector>

#include <ISO_Fortran_binding.h>

#include <gtest/gtest.h>
#include <micm/solver/chapman_ode_solver.hpp>

extern "C" {
  void Finit(void);
  void get_names( CFI_cdesc_t * );
}

TEST(RegressionChapmanODESolver, ReactionNames){
  Finit();

  micm::ChapmanODESolver solver{};

  std::vector<std::string_view> vs;
  CFI_CDESC_T(1) names;

  get_names((CFI_cdesc_t *)&names);

  for (int i = 0; i < names.dim[0].extent; i++) {
    vs.push_back(std::string_view((char *)names.base_addr).substr(i * names.elem_len, names.elem_len));
  }


}