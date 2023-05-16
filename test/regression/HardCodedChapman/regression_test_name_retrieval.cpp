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

class RegressionChapmanODESolver : public ::testing::Test
{
protected:
     virtual void SetUp()
     {      
      // runs once for each test
     }

     virtual void TearDown()
     {
     }
public:
    static void SetUpTestSuite() {
      // runs once for the entire file
      Finit();
    }

    static void TearDownTestSuite() {
    }
};

std::vector<std::string_view> extract_names(CFI_cdesc_t* names){
  std::vector<std::string_view> vs;

  for (int i = 0; i < names->dim[0].extent; i++) {
    // fortran returns the white-space padded string
    // find the first whitespace character at each address to determine the actual
    // length of the string
    // only works since whitespace characters aren't important for reaction names
    char* addr = (char *)(names->base_addr) + i * names->elem_len;
    char* first_space = strchr(addr, ' ');
    size_t strlen = first_space - addr;
    vs.push_back(std::string_view(addr).substr(0, strlen));
  }

  return vs;
}

TEST_F(RegressionChapmanODESolver, ReactionNames){
  micm::ChapmanODESolver solver{};

  CFI_CDESC_T(1) names;

  CFI_establish((CFI_cdesc_t *)&names, NULL,
                      CFI_attribute_pointer,
                      CFI_type_char, 0, (CFI_rank_t)1, NULL);

  get_reaction_names((CFI_cdesc_t *)&names);

  std::vector<std::string_view> vs = extract_names((CFI_cdesc_t *)&names);

  for(const auto& elem : solver.reaction_names())
  {
    auto it = std::find(vs.begin(), vs.end(), elem);
    EXPECT_TRUE(it != vs.end());
  }
}

TEST_F(RegressionChapmanODESolver, PhotolysisNames){
  micm::ChapmanODESolver solver{};

  CFI_CDESC_T(1) names;
  CFI_establish((CFI_cdesc_t *)&names, NULL,
                      CFI_attribute_pointer,
                      CFI_type_char, 0, (CFI_rank_t)1, NULL);

  get_photolysis_names((CFI_cdesc_t *)&names);

  std::vector<std::string_view> vs = extract_names((CFI_cdesc_t *)&names);

  for(const auto& elem : solver.photolysis_names())
  {
    auto it = std::find(vs.begin(), vs.end(), elem);
    EXPECT_TRUE(it != vs.end());
  }
}

TEST_F(RegressionChapmanODESolver, SpeciesNames){
  micm::ChapmanODESolver solver{};

  CFI_CDESC_T(1) names;
  CFI_establish((CFI_cdesc_t *)&names, NULL,
                      CFI_attribute_pointer,
                      CFI_type_char, 0, (CFI_rank_t)1, NULL);

  get_species_names((CFI_cdesc_t *)&names);

  std::vector<std::string_view> vs = extract_names((CFI_cdesc_t *)&names);

  for(const auto& elem : solver.species_names())
  {
    auto it = std::find(vs.begin(), vs.end(), elem);
    EXPECT_TRUE(it != vs.end());
  }
}