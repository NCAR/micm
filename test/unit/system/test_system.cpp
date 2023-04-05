#include <micm/system/system.hpp>
#include <micm/system/species.hpp>

#include <gtest/gtest.h>

micm::SystemParameters FullSetOfParameters(){
  micm::SystemParameters parameters;

  parameters.gas_phase_ = micm::Phase();

  parameters.phases_ = std::vector<micm::Phase> {
    micm::Phase(),
    micm::Phase(),
    micm::Phase(),
  };

  return parameters;
}

TEST(System, DefaultConstructor){
  micm::System system{};
}

TEST(System, ConstructorWithPhase){
  micm::SystemParameters parameters;
  parameters.gas_phase_ = micm::Phase();

  micm::System system{parameters};

  EXPECT_EQ(system.gas_phase_.species_.size(), 0);
}

TEST(System, ConstructorWithMultiplePhases){
  micm::SystemParameters parameters;
  parameters.phases_ = std::vector<micm::Phase> {
    micm::Phase(),
    micm::Phase(),
    micm::Phase(),
  };

  micm::System system{parameters};

  EXPECT_EQ(system.phases_.size(), 3);
}

TEST(System, ConstructorWithAllParameters){
  micm::System system{FullSetOfParameters()};

  EXPECT_EQ(system.gas_phase_.species_.size(), 0);
  EXPECT_EQ(system.phases_.size(), 3);
}