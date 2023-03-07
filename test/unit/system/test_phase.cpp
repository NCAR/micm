#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <gtest/gtest.h>

TEST(Phase, DefaultConstructor){
  micm::Phase phase{};

  EXPECT_EQ(phase.species_.size(), 0);
}

TEST(Phase, ConstructorWithVector){
  micm::Phase phase(std::vector<micm::Species>({micm::Species("species1"), micm::Species("species2")}));

  EXPECT_EQ(phase.species_.size(), 2);
}

TEST(Phase, CopyConstructor){
  micm::Phase phase(std::vector<micm::Species>({micm::Species("species1"), micm::Species("species2")}));
  micm::Phase phase2(phase);

  EXPECT_EQ(phase.species_.size(), 2);
}

TEST(Phase, CopyAssignment){
  micm::Phase phase(std::vector<micm::Species>({micm::Species("species1"), micm::Species("species2")}));
  micm::Phase phase2 = phase;

  EXPECT_EQ(phase.species_.size(), 2);
}
