#include <gtest/gtest.h>

#include <algorithm>
#include <micm/system/species.hpp>
#include <micm/system/system.hpp>

TEST(System, ConstructorWithAllParameters)
{
  std::vector<micm::Species> speciesA = { micm::Species("species1"), micm::Species("species2") };
  std::vector<micm::Species> speciesB = { micm::Species("species3"), micm::Species("species4") };

  micm::Phase phase = speciesA;
  std::unordered_map<std::string, micm::Phase> phases = { { "phase1", speciesA }, { "phase2", speciesB } };

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

  std::vector<int> reorder{ 3, 2, 1, 0, 5, 4 };
  auto reordered_names = system.UniqueNames([&](const std::vector<std::string> variables, const std::size_t i)
                                            { return variables[reorder[i]]; });
  EXPECT_EQ(reordered_names.size(), 6);
  EXPECT_EQ(reordered_names[0], names[3]);
  EXPECT_EQ(reordered_names[1], names[2]);
  EXPECT_EQ(reordered_names[2], names[1]);
  EXPECT_EQ(reordered_names[3], names[0]);
  EXPECT_EQ(reordered_names[4], names[5]);
  EXPECT_EQ(reordered_names[5], names[4]);
}

TEST(System, ConstructorWithParameterizedSpecies)
{
  std::vector<micm::Species> speciesA = { micm::Species("species1"), micm::Species("species2") };
  std::vector<micm::Species> speciesB = { micm::Species("species3"), micm::Species("species4") };
  auto param_species = micm::Species("paramSpecies");
  param_species.parameterize_ = [](const micm::Conditions& c) { return 64.2; };
  speciesA.push_back(param_species);

  micm::Phase phase = speciesA;
  std::unordered_map<std::string, micm::Phase> phases = { { "phase1", speciesA }, { "phase2", speciesB } };

  micm::System system = { micm::SystemParameters{ .gas_phase_ = phase, .phases_ = phases } };

  EXPECT_EQ(system.gas_phase_.species_.size(), 3);
  EXPECT_EQ(system.phases_.size(), 2);
  EXPECT_EQ(system.StateSize(), 6);

  auto names = system.UniqueNames();

  EXPECT_EQ(names.size(), 6);
  EXPECT_NE(std::find(names.begin(), names.end(), "species1"), names.end());
  EXPECT_NE(std::find(names.begin(), names.end(), "species2"), names.end());
  EXPECT_NE(std::find(names.begin(), names.end(), "phase1.species1"), names.end());
  EXPECT_NE(std::find(names.begin(), names.end(), "phase1.species2"), names.end());
  EXPECT_NE(std::find(names.begin(), names.end(), "phase2.species3"), names.end());
  EXPECT_NE(std::find(names.begin(), names.end(), "phase2.species4"), names.end());

  std::vector<int> reorder{ 3, 2, 1, 0, 5, 4 };
  auto reordered_names = system.UniqueNames([&](const std::vector<std::string> variables, const std::size_t i)
                                            { return variables[reorder[i]]; });
  EXPECT_EQ(reordered_names.size(), 6);
  EXPECT_EQ(reordered_names[0], names[3]);
  EXPECT_EQ(reordered_names[1], names[2]);
  EXPECT_EQ(reordered_names[2], names[1]);
  EXPECT_EQ(reordered_names[3], names[0]);
  EXPECT_EQ(reordered_names[4], names[5]);
  EXPECT_EQ(reordered_names[5], names[4]);
}