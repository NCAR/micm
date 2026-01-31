#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

#include <gtest/gtest.h>

using namespace micm;

TEST(PhaseSpecies, ConstructsWithoutDiffusionCoefficient)
{
  Species CO2("CO2");
  PhaseSpecies gas_CO2(CO2);
  EXPECT_EQ(gas_CO2.species_.name_, "CO2");
  EXPECT_FALSE(gas_CO2.diffusion_coefficient_.has_value());
}

TEST(PhaseSpecies, ConstructsWithDiffusionCoefficient)
{
  Species CO2("CO2");
  double diff_coeff = 1.23e-5;
  PhaseSpecies gas_CO2(CO2, diff_coeff);
  EXPECT_EQ(gas_CO2.species_.name_, "CO2");
  ASSERT_TRUE(gas_CO2.diffusion_coefficient_.has_value());
  EXPECT_DOUBLE_EQ(gas_CO2.diffusion_coefficient_.value(), diff_coeff);
}

TEST(PhaseSpecies, SetDiffusionCoefficient)
{
  Species CO2("CO2");
  PhaseSpecies gas_CO2(CO2);
  EXPECT_FALSE(gas_CO2.diffusion_coefficient_.has_value());

  double diff_coeff = 2.5e-6;
  gas_CO2.SetDiffusionCoefficient(diff_coeff);
  ASSERT_TRUE(gas_CO2.diffusion_coefficient_.has_value());
  EXPECT_DOUBLE_EQ(gas_CO2.diffusion_coefficient_.value(), diff_coeff);
}

TEST(Phase, Constructor)
{
  Phase phase(
      "test_phase", std::vector<PhaseSpecies>({ PhaseSpecies(Species("species1")), PhaseSpecies(Species("species2")) }));
  EXPECT_EQ(phase.name_, "test_phase");
  EXPECT_EQ(phase.phase_species_.size(), 2);
  EXPECT_EQ(phase.StateSize(), 2);

  auto names = phase.UniqueNames();
  EXPECT_EQ(names[0], "test_phase.species1");
  EXPECT_EQ(names[1], "test_phase.species2");
}

TEST(Phase, ConstructorWithParameterizedSpecies)
{
  Species foo("foo");
  Species bar("bar");
  Species baz("baz");
  PhaseSpecies gas_foo(foo);
  PhaseSpecies gas_bar(bar);
  PhaseSpecies gas_baz(baz);

  gas_bar.species_.parameterize_ = [](const Conditions& c) { return 42.0; };
  Phase phase("gas", std::vector<PhaseSpecies>({ gas_foo, gas_bar, gas_baz }));

  EXPECT_EQ(phase.phase_species_.size(), 3);
  EXPECT_EQ(phase.StateSize(), 2);
}

TEST(Phase, UniqueNamesWithParameterizedSpecies)
{
  Species foo("foo");
  Species bar("bar");
  Species baz("baz");
  PhaseSpecies gas_foo(foo);
  PhaseSpecies gas_bar(bar);
  PhaseSpecies gas_baz(baz);

  gas_bar.species_.parameterize_ = [](const Conditions& c) { return 42.0; };
  Phase phase("gas", std::vector<PhaseSpecies>({ gas_foo, gas_bar, gas_baz }));

  auto names = phase.UniqueNames();
  EXPECT_EQ(names.size(), 2);
  EXPECT_EQ(names[0], "gas.foo");
  EXPECT_EQ(names[1], "gas.baz");
}