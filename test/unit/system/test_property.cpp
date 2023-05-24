#include <gtest/gtest.h>

#include <micm/system/property.hpp>

TEST(Property, Constructor)
{
  micm::Property property("name", "units", 1.0);

  EXPECT_EQ(property.name_, "name");
  EXPECT_EQ(property.units_, "units");
  EXPECT_EQ(property.value_, 1.0);
}

TEST(Property, CopyConstructor)
{
  micm::Property property("name", "units", 1.0);
  micm::Property property2(property);

  EXPECT_EQ(property2.name_, "name");
  EXPECT_EQ(property2.units_, "units");
  EXPECT_EQ(property2.value_, 1.0);
}

TEST(Property, CopyAssignment)
{
  micm::Property property("name", "units", 1.0);
  micm::Property property2 = property;

  EXPECT_EQ(property2.name_, "name");
  EXPECT_EQ(property2.units_, "units");
  EXPECT_EQ(property2.value_, 1.0);
}
