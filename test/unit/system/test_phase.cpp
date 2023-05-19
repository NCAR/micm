#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <gtest/gtest.h>

TEST(Phase, ConstructorWithVector){
  micm::Phase phase(std::vector<micm::Species>({micm::Species("species1"), micm::Species("species2")}));

  EXPECT_EQ(phase.species_.size(), 2);
}