#include <string_view>
#include <vector>
#include <algorithm>

#include <ISO_Fortran_binding.h>

#include <gtest/gtest.h>
#include <micm/solver/chapman_ode_solver.hpp>

extern "C" {
  void Finit(void);
  void get_reaction_names( CFI_cdesc_t * );
  void get_photolysis_names( CFI_cdesc_t * );
  void get_species_names( CFI_cdesc_t * );
}

std::vector<std::string_view> extract_names(CFI_cdesc_t* names){
  std::vector<std::string_view> vs;

  for (int i = 0; i < names->dim[0].extent; i++) {
    // fortran returns the white-space padded string
    // find the first whitespace character at each address to determine the actual
    // length of the string
    // only works since whitespace characters aren't important for reaction names
    char* addr = (char *)(names->base_addr + i * names->elem_len);
    char* first_space = strchr(addr, ' ');
    size_t strlen = first_space - addr;
    vs.push_back(std::string_view(addr).substr(0, strlen));
  }

  return vs;
}

TEST(RegressionChapmanODESolver, ReactionNames){
  Finit();
  micm::ChapmanODESolver solver{};

  CFI_CDESC_T(1) names;

  get_reaction_names((CFI_cdesc_t *)&names);

  std::vector<std::string_view> vs = extract_names((CFI_cdesc_t *)&names);

  for(const auto& elem : solver.reaction_names())
  {
    auto it = std::find(vs.begin(), vs.end(), elem);
    EXPECT_TRUE(it != vs.end());
  }
}

// TEST(RegressionChapmanODESolver, PhotolysisNames){
//   Finit();
//   micm::ChapmanODESolver solver{};

//   CFI_CDESC_T(1) names;

//   get_photolysis_names((CFI_cdesc_t *)&names);

//   std::vector<std::string_view> vs = extract_names((CFI_cdesc_t *)&names);

//   for(const auto& elem : solver.photolysis_names())
//   {
//     auto it = std::find(vs.begin(), vs.end(), elem);
//     EXPECT_TRUE(it != vs.end());
//   }
// }

// TEST(RegressionChapmanODESolver, SpeciesNames){
//   Finit();
//   micm::ChapmanODESolver solver{};

//   CFI_CDESC_T(1) names;

//   get_species_names((CFI_cdesc_t *)&names);

//   std::vector<std::string_view> vs = extract_names((CFI_cdesc_t *)&names);

//   for(const auto& elem : solver.species_names())
//   {
//     auto it = std::find(vs.begin(), vs.end(), elem);
//     EXPECT_TRUE(it != vs.end());
//   }
// }