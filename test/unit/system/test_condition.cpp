#include <micm/system/condition.hpp>

#include <gtest/gtest.h>

TEST(Condition, ConstructorWithArguments){
  micm::Condition condition("name", "units");

  EXPECT_EQ(condition.name_, "name");
  EXPECT_EQ(condition.units_, "units");
}

TEST(Condition, CopyAssignment){
  micm::Condition condition("name", "units");
  micm::Condition condition2 = condition;

  EXPECT_EQ(condition2.name_, "name");
  EXPECT_EQ(condition2.units_, "units");
}

TEST(Condition, CopyConstructor){
  micm::Condition condition("name", "units");
  micm::Condition condition2(condition);

  EXPECT_EQ(condition2.name_, "name");
  EXPECT_EQ(condition2.units_, "units");
}