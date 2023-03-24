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

  parameters.conditions_ = { micm::Condition("name", "units") };

  parameters.photolysis_reactions_ = micm::IntraPhaseProcess<micm::PhotolysisRateConstant>(
    std::vector<micm::Species> { micm::Species("O2") },
    std::vector<micm::Species> { micm::Species("O"), micm::Species("O") },
    micm::PhotolysisRateConstant(1.1)
  );

  micm::ArrheniusRateConstantParameters arrhenius_parameters;
  arrhenius_parameters.A_ = 3.3e-11;
  arrhenius_parameters.C_ = 55;
  parameters.arrhenius_reactions_ = micm::IntraPhaseProcess<micm::ArrheniusRateConstant>(
    std::vector<micm::Species> { micm::Species("O1D"), micm::Species("O2") },
    std::vector<micm::Species> { micm::Species("O"), micm::Species("O2") },
    micm::ArrheniusRateConstant(arrhenius_parameters)
  );

  return parameters;
}

void TestFullSetOfParameters(micm::System system){
  EXPECT_EQ(system.gas_phase_.species_.size(), 0);
  EXPECT_EQ(system.phases_.size(), 3);
  EXPECT_EQ(system.conditions_.size(), 1);
  EXPECT_EQ(system.conditions_[0].name_, "name");
  EXPECT_EQ(system.conditions_[0].units_, "units");
  EXPECT_EQ(system.photolysis_reactions_.rate_.rate_, 1.1);
  EXPECT_EQ(system.photolysis_reactions_.reactants_.size(), 1);
  EXPECT_EQ(system.photolysis_reactions_.products_.size(), 2);
  EXPECT_EQ(system.arrhenius_reactions_.rate_.parameters_.A_, 3.3e-11);
  EXPECT_EQ(system.arrhenius_reactions_.rate_.parameters_.C_, 55);
  EXPECT_EQ(system.arrhenius_reactions_.reactants_.size(), 2);
  EXPECT_EQ(system.arrhenius_reactions_.products_.size(), 2);
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

TEST(System, ConstructorWithCondition){
  micm::SystemParameters parameters;
  parameters.conditions_ = { micm::Condition("name", "units") };

  micm::System system(parameters);

  EXPECT_EQ(system.conditions_.size(), 1);
  EXPECT_EQ(system.conditions_[0].name_, "name");
  EXPECT_EQ(system.conditions_[0].units_, "units");
}

TEST(System, ConstructorWithPhotolysisReactions){
  micm::SystemParameters parameters;
  parameters.photolysis_reactions_ = micm::IntraPhaseProcess<micm::PhotolysisRateConstant>(
    std::vector<micm::Species> { micm::Species("O2") },
    std::vector<micm::Species> { micm::Species("O"), micm::Species("O") },
    micm::PhotolysisRateConstant(1.1)
  );

  micm::System system(parameters);

  EXPECT_EQ(system.photolysis_reactions_.rate_.rate_, 1.1);
  EXPECT_EQ(system.photolysis_reactions_.reactants_.size(), 1);
  EXPECT_EQ(system.photolysis_reactions_.products_.size(), 2);
}

TEST(System, ConstructorWithArrheniusReaction){
  micm::SystemParameters parameters;
  micm::ArrheniusRateConstantParameters arrhenius_parameters;
  arrhenius_parameters.A_ = 3.3e-11;
  arrhenius_parameters.C_ = 55;

  parameters.arrhenius_reactions_ = micm::IntraPhaseProcess<micm::ArrheniusRateConstant>(
    std::vector<micm::Species> { micm::Species("O1D"), micm::Species("O2") },
    std::vector<micm::Species> { micm::Species("O"), micm::Species("O2") },
    micm::ArrheniusRateConstant(arrhenius_parameters)
  );

  micm::System system(parameters);

  EXPECT_EQ(system.arrhenius_reactions_.rate_.parameters_.A_, 3.3e-11);
  EXPECT_EQ(system.arrhenius_reactions_.rate_.parameters_.C_, 55);
  EXPECT_EQ(system.arrhenius_reactions_.reactants_.size(), 2);
  EXPECT_EQ(system.arrhenius_reactions_.products_.size(), 2);
}

TEST(System, CopyConstructor){
  micm::System system{FullSetOfParameters()};
  auto system2{system};

  TestFullSetOfParameters(system2);
}

TEST(System, MoveConstructor){
  micm::System system{micm::System(FullSetOfParameters())};
  TestFullSetOfParameters(system);
}
