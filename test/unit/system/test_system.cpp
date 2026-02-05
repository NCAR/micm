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
  PhaseSpecies gas_foo(foo);
  PhaseSpecies gas_bar(bar);

  Phase gas_phase("gas", std::vector<PhaseSpecies>({ gas_foo, gas_bar }));
  System system = { SystemParameters{ .gas_phase_ = gas_phase } };

  EXPECT_EQ(system.gas_phase_.phase_species_.size(), 2);

  auto names = system.UniqueNames();
  std::vector<std::string> expected = { "foo", "bar" };
  std::multiset<std::string> name_set(names.begin(), names.end());
  std::multiset<std::string> expected_set(expected.begin(), expected.end());

  EXPECT_EQ(names.size(), expected.size());
  EXPECT_EQ(name_set, expected_set);
}

TEST(System, ConstructorWithParameterizedSpecies)
{
  Species foo("foo");
  Species bar("bar");
  PhaseSpecies gas_foo(foo);
  PhaseSpecies gas_bar(bar);

  // Parameterized species
  Species param_species = Species("paramSpecies");
  param_species.parameterize_ = [](const Conditions& c) { return 64.2; };
  PhaseSpecies gas_param_species(param_species);

  Phase gas_phase("gas", std::vector<PhaseSpecies>({ gas_foo, gas_bar, gas_param_species }));
  System system = { SystemParameters{ .gas_phase_ = gas_phase } };

  EXPECT_EQ(system.gas_phase_.phase_species_.size(), 3);
  EXPECT_EQ(system.StateSize(), 3 - 1);  // One parameterized species

  auto names = system.UniqueNames();
  std::vector<std::string> expected = { "foo", "bar" };
  std::multiset<std::string> name_set(names.begin(), names.end());
  std::multiset<std::string> expected_set(expected.begin(), expected.end());

  EXPECT_EQ(names.size(), expected.size());
  EXPECT_EQ(name_set, expected_set);

  std::vector<int> reorder{ 1, 0 };
  auto reordered_names = system.UniqueNames([&](const std::vector<std::string> variables, const std::size_t i)
                                            { return variables[reorder[i]]; });
  EXPECT_EQ(reordered_names.size(), 2);
  EXPECT_EQ(reordered_names[0], names[1]);
  EXPECT_EQ(reordered_names[1], names[0]);
}
