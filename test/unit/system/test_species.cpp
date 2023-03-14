#include <micm/system/species.hpp>

#include <gtest/gtest.h>

TEST(Species, DefaultConstructor){
  micm::Species species{};

  EXPECT_EQ(species.properties_.size(), 0);
}

TEST(Species, StringConstructor){
  micm::Species species("thing");

  EXPECT_EQ(species.name_, "thing");
}

TEST(Species, StringAndVectorConstructor){
  micm::Species species("thing", std::vector<micm::Property>{
    micm::Property("name", "units", 1.0),
    micm::Property("name2", "units2", 2.0)
  });

  EXPECT_EQ(species.name_, "thing");
  EXPECT_EQ(species.properties_.size(), 2);
  EXPECT_EQ(species.properties_[0].name_, "name");
  EXPECT_EQ(species.properties_[0].units_, "units");
  EXPECT_EQ(species.properties_[0].value_, 1.0);

  EXPECT_EQ(species.properties_[1].name_, "name2");
  EXPECT_EQ(species.properties_[1].units_, "units2");
  EXPECT_EQ(species.properties_[1].value_, 2.0);
}

TEST(Species, CopyConstructor){
  micm::Species species("thing", std::vector<micm::Property>{
    micm::Property("name", "units", 1.0),
    micm::Property("name2", "units2", 2.0)
  });

  micm::Species species2(species);

  EXPECT_EQ(species2.name_, "thing");
  EXPECT_EQ(species2.properties_.size(), 2);
  EXPECT_EQ(species2.properties_[0].name_, "name");
  EXPECT_EQ(species2.properties_[0].units_, "units");
  EXPECT_EQ(species2.properties_[0].value_, 1.0);

  EXPECT_EQ(species2.properties_[1].name_, "name2");
  EXPECT_EQ(species2.properties_[1].units_, "units2");
  EXPECT_EQ(species2.properties_[1].value_, 2.0);
}

TEST(Species, CopyAssignment){
  micm::Species species("thing", std::vector<micm::Property>{
    micm::Property("name", "units", 1.0),
    micm::Property("name2", "units2", 2.0)
  });

  micm::Species species2 = species;

  EXPECT_EQ(species2.name_, "thing");
  EXPECT_EQ(species2.properties_.size(), 2);
  EXPECT_EQ(species2.properties_[0].name_, "name");
  EXPECT_EQ(species2.properties_[0].units_, "units");
  EXPECT_EQ(species2.properties_[0].value_, 1.0);

  EXPECT_EQ(species2.properties_[1].name_, "name2");
  EXPECT_EQ(species2.properties_[1].units_, "units2");
  EXPECT_EQ(species2.properties_[1].value_, 2.0);
}
