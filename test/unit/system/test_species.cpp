#include <gtest/gtest.h>

#include <micm/system/species.hpp>

TEST(Species, StringConstructor)
{
  micm::Species species("thing");
  EXPECT_EQ(species.name_, "thing");
  EXPECT_FALSE(species.IsParameterized());
  species.parameterize_ = [](const micm::Conditions& c) { return c.temperature_ * 100.0; };
  EXPECT_TRUE(species.IsParameterized());
  EXPECT_EQ(species.parameterize_({ .temperature_ = 12.5 }), 12.5 * 100.0);
}

TEST(Species, StringAndVectorConstructor)
{
  micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });

  EXPECT_EQ(species.name_, "thing");
  EXPECT_EQ(species.properties_.size(), 2);
  EXPECT_EQ(species.properties_["name [units]"], 1.0);
  EXPECT_EQ(species.properties_["name2 [units2]"], 2.0);
}

TEST(Species, ThirdBody)
{
  micm::Species species = micm::Species::ThirdBody();
  EXPECT_EQ(species.name_, "M");
  EXPECT_TRUE(species.IsParameterized());
  EXPECT_EQ(species.parameterize_({ .air_density_ = 42.4 }), 42.4);
}

TEST(Species, CopyConstructor)
{
  {
    micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });

    micm::Species species2(species);

    EXPECT_EQ(species2.name_, "thing");
    EXPECT_EQ(species2.properties_.size(), 2);
    EXPECT_EQ(species2.properties_["name [units]"], 1.0);
    EXPECT_EQ(species2.properties_["name2 [units2]"], 2.0);
    EXPECT_FALSE(species2.IsParameterized());
  }
  {
    micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });
    species.parameterize_ = [](const micm::Conditions& c) { return 15.4; };

    micm::Species species2(species);

    EXPECT_EQ(species2.name_, "thing");
    EXPECT_EQ(species2.properties_.size(), 2);
    EXPECT_EQ(species2.properties_["name [units]"], 1.0);
    EXPECT_EQ(species2.properties_["name2 [units2]"], 2.0);
    EXPECT_TRUE(species2.IsParameterized());
    EXPECT_EQ(species.parameterize_({}), 15.4);
  }
}

TEST(Species, CopyAssignment)
{
  {
    micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });

    micm::Species species2 = species;

    EXPECT_EQ(species2.name_, "thing");
    EXPECT_EQ(species2.properties_.size(), 2);
    EXPECT_EQ(species2.properties_["name [units]"], 1.0);
    EXPECT_EQ(species2.properties_["name2 [units2]"], 2.0);
    EXPECT_FALSE(species2.IsParameterized());
  }
  {
    micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });
    species.parameterize_ = [](const micm::Conditions& c) { return 15.4; };

    micm::Species species2 = species;

    EXPECT_EQ(species2.name_, "thing");
    EXPECT_EQ(species2.properties_.size(), 2);
    EXPECT_EQ(species2.properties_["name [units]"], 1.0);
    EXPECT_EQ(species2.properties_["name2 [units2]"], 2.0);
    EXPECT_TRUE(species2.IsParameterized());
    EXPECT_EQ(species.parameterize_({}), 15.4);
  }
}
