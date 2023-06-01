#include <gtest/gtest.h>

#include <micm/system/species.hpp>
#include <micm/system/system.hpp>

TEST(System, ConstructorWithAllParameters)
{
  std::vector<micm::Species> speciesA = { micm::Species("species1"), micm::Species("species2") };
  std::vector<micm::Species> speciesB = { micm::Species("species3"), micm::Species("species4") };

  micm::Phase phase = speciesA;
  std::unordered_map<std::string, micm::Phase> phases = {
    { "phase1", speciesA },
    { "phase2", speciesB }
  };

  micm::System system = { micm::SystemParameters{ .gas_phase_ = phase, .phases_ = phases } };

  EXPECT_EQ(system.gas_phase_.species_.size(), 2);
  EXPECT_EQ(system.phases_.size(), 2);

  auto names = system.UniqueNames();

  EXPECT_EQ(names.size(), 6);
  EXPECT_NE(std::find(names.begin(), names.end(), "species1"), names.end());
  EXPECT_NE(std::find(names.begin(), names.end(), "species2"), names.end());
  EXPECT_NE(std::find(names.begin(), names.end(), "phase1.species1"), names.end());
  EXPECT_NE(std::find(names.begin(), names.end(), "phase1.species2"), names.end());
  EXPECT_NE(std::find(names.begin(), names.end(), "phase2.species3"), names.end());
  EXPECT_NE(std::find(names.begin(), names.end(), "phase2.species4"), names.end());
}