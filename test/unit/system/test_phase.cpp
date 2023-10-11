#include <gtest/gtest.h>

#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>

TEST(Phase, ConstructorWithVector)
{
  micm::Phase phase(std::vector<micm::Species>({ micm::Species("species1"), micm::Species("species2") }));

  EXPECT_EQ(phase.species_.size(), 2);
  EXPECT_EQ(phase.StateSize(), 2);

  auto names = phase.UniqueNames();
  EXPECT_EQ(names[0], "species1");
  EXPECT_EQ(names[1], "species2");
}

TEST(Phase, ConstructorWithParameterizedSpecies)
{
  auto foo = micm::Species("foo");
  auto bar = micm::Species("bar");
  auto baz = micm::Species("baz");

  bar.parameterize_ = [](const micm::Conditions& c) { return 42.0; };
  micm::Phase phase(std::vector<micm::Species>({ foo, bar, baz }));

  EXPECT_EQ(phase.species_.size(), 3);
  EXPECT_EQ(phase.StateSize(), 2);

  auto names = phase.UniqueNames();
  EXPECT_EQ(names[0], "foo");
  EXPECT_EQ(names[1], "baz");
}