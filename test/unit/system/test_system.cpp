#include <micm/system/system.hpp>

#include <gtest/gtest.h>

TEST(System, DefaultConstructor){
  micm::System system{};
}

TEST(System, ConstructorWithPhase){
  micm::System system{micm::Phase()};
  EXPECT_EQ(system.gas_phase_.species_.size(), 0);
}

TEST(System, ConstructorWithPhaseAndCondition){
  micm::System system(micm::Phase(), micm::Condition("name", "units"));
  EXPECT_EQ(system.gas_phase_.species_.size(), 0);
  EXPECT_EQ(system.conditions_.size(), 1);
  EXPECT_EQ(system.conditions_[0].name_, "name");
  EXPECT_EQ(system.conditions_[0].units_, "units");
}

TEST(System, ConstructorWithPhaseAndVectorCondition){
  micm::System system(micm::Phase(), std::vector<micm::Condition>{
    micm::Condition("name", "units"),
    micm::Condition("name2", "units2")
  });
  EXPECT_EQ(system.gas_phase_.species_.size(), 0);
  EXPECT_EQ(system.conditions_.size(), 2);
  EXPECT_EQ(system.conditions_[0].name_, "name");
  EXPECT_EQ(system.conditions_[0].units_, "units");
  EXPECT_EQ(system.conditions_[1].name_, "name2");
  EXPECT_EQ(system.conditions_[1].units_, "units2");
}
