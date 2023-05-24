#include <gtest/gtest.h>

#include <micm/system/species.hpp>
#include <micm/system/system.hpp>

TEST(System, ConstructorWithAllParameters)
{
  std::vector<micm::Species> speciesA = { micm::Species("species1"), micm::Species("species2") };
  std::vector<micm::Species> speciesB = { micm::Species("species3"), micm::Species("species4") };

  micm::Phase phase = speciesA;
  std::vector<micm::Phase> phases = { speciesA, speciesB };

  micm::System system = { micm::SystemParameters{ .gas_phase_ = phase, .phases_ = phases } };

  EXPECT_EQ(system.gas_phase_.species_.size(), 2);
  EXPECT_EQ(system.phases_.size(), 2);
}