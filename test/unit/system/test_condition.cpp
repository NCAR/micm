#include <micm/system/condition.hpp>

#include <gtest/gtest.h>

TEST(Condition, DefaultConstructor){
  micm::Condition condition{};
}

TEST(Condition, ConstructorWithArguments){
  micm::Condition condition("name", "units");
}