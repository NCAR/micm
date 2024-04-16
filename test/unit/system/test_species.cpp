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
  EXPECT_EQ(species.properties_double_.size(), 2);
  EXPECT_EQ(species.GetProperty<double>("name [units]"), 1.0);
  EXPECT_EQ(species.GetProperty<double>("name2 [units2]"), 2.0);
}

TEST(Species, GetProperty)
{
  micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });

  EXPECT_EQ(species.GetProperty<double>("name [units]"), 1.0);
  EXPECT_EQ(species.GetProperty<double>("name2 [units2]"), 2.0);
  EXPECT_THROW(
      {
        try
        {
          species.GetProperty<std::string>("not there");
        }
        catch (std::runtime_error& e)
        {
          EXPECT_STREQ(e.what(), "Species property 'not there' not found");
          throw;
        }
      },
      std::runtime_error);
  EXPECT_THROW(
      {
        try
        {
          species.GetProperty<double>("not there");
        }
        catch (std::runtime_error& e)
        {
          EXPECT_STREQ(e.what(), "Species property 'not there' not found");
          throw;
        }
      },
      std::runtime_error);
  EXPECT_THROW(
      {
        try
        {
          species.GetProperty<int>("not there");
        }
        catch (std::runtime_error& e)
        {
          EXPECT_STREQ(e.what(), "Species property 'not there' not found");
          throw;
        }
      },
      std::runtime_error);
  EXPECT_THROW(
      {
        try
        {
          species.GetProperty<bool>("not there");
        }
        catch (std::runtime_error& e)
        {
          EXPECT_STREQ(e.what(), "Species property 'not there' not found");
          throw;
        }
      },
      std::runtime_error);
  EXPECT_THROW(
      {
        try
        {
          species.GetProperty<long double>("name [units]");
        }
        catch (std::runtime_error& e)
        {
          EXPECT_STREQ(e.what(), "Invalid type for species property");
          throw;
        }
      },
      std::runtime_error);
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
    species.SetProperty("foo", "bar");
    species.SetProperty("baz", 42);
    species.SetProperty("qux", true);

    micm::Species species2(species);

    EXPECT_EQ(species2.name_, "thing");
    EXPECT_EQ(species2.properties_double_.size(), 2);
    EXPECT_EQ(species2.GetProperty<double>("name [units]"), 1.0);
    EXPECT_EQ(species2.GetProperty<double>("name2 [units2]"), 2.0);
    EXPECT_EQ(species2.GetProperty<std::string>("foo"), "bar");
    EXPECT_EQ(species2.GetProperty<int>("baz"), 42);
    EXPECT_EQ(species2.GetProperty<bool>("qux"), true);
    EXPECT_FALSE(species2.IsParameterized());
  }
  {
    micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });
    species.parameterize_ = [](const micm::Conditions& c) { return 15.4; };
    species.SetProperty("foo", "bar");
    species.SetProperty("baz", 42);
    species.SetProperty("qux", true);

    micm::Species species2(species);

    EXPECT_EQ(species2.name_, "thing");
    EXPECT_EQ(species2.properties_double_.size(), 2);
    EXPECT_EQ(species2.GetProperty<double>("name [units]"), 1.0);
    EXPECT_EQ(species2.GetProperty<double>("name2 [units2]"), 2.0);
    EXPECT_EQ(species2.GetProperty<std::string>("foo"), "bar");
    EXPECT_EQ(species2.GetProperty<int>("baz"), 42);
    EXPECT_EQ(species2.GetProperty<bool>("qux"), true);
    EXPECT_TRUE(species2.IsParameterized());
    EXPECT_EQ(species.parameterize_({}), 15.4);
  }
}

TEST(Species, CopyAssignment)
{
  {
    micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });
    species.SetProperty("foo", "bar");
    species.SetProperty("baz", 42);
    species.SetProperty("qux", true);

    micm::Species species2 = species;

    EXPECT_EQ(species2.name_, "thing");
    EXPECT_EQ(species2.properties_double_.size(), 2);
    EXPECT_EQ(species2.GetProperty<double>("name [units]"), 1.0);
    EXPECT_EQ(species2.GetProperty<double>("name2 [units2]"), 2.0);
    EXPECT_EQ(species2.GetProperty<std::string>("foo"), "bar");
    EXPECT_EQ(species2.GetProperty<int>("baz"), 42);
    EXPECT_EQ(species2.GetProperty<bool>("qux"), true);
    EXPECT_FALSE(species2.IsParameterized());
  }
  {
    micm::Species species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } });
    species.parameterize_ = [](const micm::Conditions& c) { return 15.4; };
    species.SetProperty("foo", "bar");
    species.SetProperty("baz", 42);
    species.SetProperty("qux", true);

    micm::Species species2 = species;

    EXPECT_EQ(species2.name_, "thing");
    EXPECT_EQ(species2.properties_double_.size(), 2);
    EXPECT_EQ(species2.GetProperty<double>("name [units]"), 1.0);
    EXPECT_EQ(species2.GetProperty<double>("name2 [units2]"), 2.0);
    EXPECT_EQ(species2.GetProperty<std::string>("foo"), "bar");
    EXPECT_EQ(species2.GetProperty<int>("baz"), 42);
    EXPECT_EQ(species2.GetProperty<bool>("qux"), true);
    EXPECT_TRUE(species2.IsParameterized());
    EXPECT_EQ(species.parameterize_({}), 15.4);
  }
}
