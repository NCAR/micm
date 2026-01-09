#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/system/system.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <set>

using namespace micm;

TEST(System, ConstructorWithAllParameters)
{
  Species foo("foo");
  Species bar("bar");
  Species baz("baz");
  Species quiz("quiz");
  PhaseSpecies gas_foo(foo);
  PhaseSpecies gas_bar(bar);
  PhaseSpecies organic_foo(foo);
  PhaseSpecies organic_bar(bar);
  PhaseSpecies aqueous_baz(baz);
  PhaseSpecies aqueous_quiz(quiz);

  Phase gas_phase("gas", std::vector<PhaseSpecies>({ gas_foo, gas_bar }));
  Phase organic_phase("organic", std::vector<PhaseSpecies>({ organic_foo, organic_bar }));
  Phase aqueous_phase("aqueous", std::vector<PhaseSpecies>({ aqueous_baz, aqueous_quiz }));
  std::unordered_map<std::string, Phase> phases = {
    { "organic", organic_phase },
    { "aqueous", aqueous_phase },
  };
  System system = { SystemParameters{ .gas_phase_ = gas_phase, .phases_ = phases } };

  EXPECT_EQ(system.gas_phase_.phase_species_.size(), 2);
  EXPECT_EQ(system.phases_.size(), 2);
  EXPECT_EQ(system.phases_["organic"].phase_species_.size(), 2);
  EXPECT_EQ(system.phases_["aqueous"].phase_species_.size(), 2);

  auto names = system.UniqueNames();
  std::vector<std::string> expected = { "foo", "bar", "organic.foo", "organic.bar", "aqueous.baz", "aqueous.quiz" };
  std::multiset<std::string> name_set(names.begin(), names.end());
  std::multiset<std::string> expected_set(expected.begin(), expected.end());

  EXPECT_EQ(names.size(), expected.size());
  EXPECT_EQ(name_set, expected_set);
}

TEST(System, ConstructorWithParameterizedSpecies)
{
  Species foo("foo");
  Species bar("bar");
  Species baz("baz");
  Species quiz("quiz");
  PhaseSpecies gas_foo(foo);
  PhaseSpecies gas_bar(bar);
  PhaseSpecies organic_foo(foo);
  PhaseSpecies organic_bar(bar);
  PhaseSpecies aqueous_baz(baz);
  PhaseSpecies aqueous_quiz(quiz);

  // Parameterized species
  Species param_species = Species("paramSpecies");
  param_species.parameterize_ = [](const Conditions& c) { return 64.2; };
  PhaseSpecies gas_param_species(param_species);

  Phase gas_phase("gas", std::vector<PhaseSpecies>({ gas_foo, gas_bar, gas_param_species }));
  Phase organic_phase("organic", std::vector<PhaseSpecies>({ organic_foo, organic_bar }));
  Phase aqueous_phase("aqueous", std::vector<PhaseSpecies>({ aqueous_baz, aqueous_quiz }));
  std::unordered_map<std::string, Phase> phases = {
    { "organic", organic_phase },
    { "aqueous", aqueous_phase },
  };
  System system = { SystemParameters{ .gas_phase_ = gas_phase, .phases_ = phases } };

  EXPECT_EQ(system.gas_phase_.phase_species_.size(), 3);
  EXPECT_EQ(system.phases_.size(), 2);
  EXPECT_EQ(system.phases_["organic"].phase_species_.size(), 2);
  EXPECT_EQ(system.phases_["aqueous"].phase_species_.size(), 2);
  EXPECT_EQ(system.StateSize(), 6);

  auto names = system.UniqueNames();
  std::vector<std::string> expected = { "foo", "bar", "organic.foo", "organic.bar", "aqueous.baz", "aqueous.quiz" };
  std::multiset<std::string> name_set(names.begin(), names.end());
  std::multiset<std::string> expected_set(expected.begin(), expected.end());

  EXPECT_EQ(names.size(), expected.size());
  EXPECT_EQ(name_set, expected_set);

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

TEST(System, OthersIsStoredAndAccessible)
{
  SystemParameters params;
  params.others_.emplace_back("modal.aitken.number_concentration");
  params.others_.emplace_back("modal.accumulation.number_concentration");
  System sys(params);
  auto names = sys.UniqueNames();
  
  EXPECT_NE(std::find(names.begin(), names.end(), "modal.aitken.number_concentration"), names.end());
  EXPECT_NE(std::find(names.begin(), names.end(), "modal.accumulation.number_concentration"), names.end());
}