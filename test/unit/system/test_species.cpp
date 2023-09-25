#include <gtest/gtest.h>

#include <micm/system/species.hpp>

TEST(Species, StringConstructor)
{
  micm::Species species("thing");

  EXPECT_EQ(species.name_, "thing");
}

TEST(Species, StringAndVectorConstructor)
{
  micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });

  EXPECT_EQ(species.name_, "thing");
  EXPECT_EQ(species.properties_.size(), 2);
  EXPECT_EQ(species.properties_["name [units]"], 1.0);
  EXPECT_EQ(species.properties_["name2 [units2]"], 2.0);
}

TEST(Species, CopyConstructor)
{
  micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });

  micm::Species species2(species);

  EXPECT_EQ(species2.name_, "thing");
  EXPECT_EQ(species2.properties_.size(), 2);
  EXPECT_EQ(species2.properties_["name [units]"], 1.0);
  EXPECT_EQ(species2.properties_["name2 [units2]"], 2.0);
}

TEST(Species, CopyAssignment)
{
  micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });

  micm::Species species2 = species;

  EXPECT_EQ(species2.name_, "thing");
  EXPECT_EQ(species2.properties_.size(), 2);
  EXPECT_EQ(species2.properties_["name [units]"], 1.0);
  EXPECT_EQ(species2.properties_["name2 [units2]"], 2.0);
}
