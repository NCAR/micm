#include <micm/system/system.hpp>
#include <micm/system/species.hpp>

#include <gtest/gtest.h>

TEST(Phase, ConstructorWithVector){
  micm::Phase phase(std::vector<micm::Species>({micm::Species("species1"), micm::Species("species2")}));

  EXPECT_EQ(phase.species_.size(), 2);
}

//  TestFullSetOfParameters
// class TestSystemParameters
// {
//   public:
//     static inline std::vector<micm::Species> speciesA = {micm::Species("species1"), micm::Species("species2")};
//     static inline std::vector<micm::Species> speciesB = {micm::Species("species3"), micm::Species("species4")};
    
//     static inline micm::Phase phase = speciesA;
//     static inline std::vector<micm::Phase> phases = {speciesA, speciesB};
    
//     static inline micm::System system = {phase, phases};
// };


// micm::SystemParameters FullSetOfParameters(){

//   const Phase gas_phase_;
//   /// @brief This is a catchall for anything that is not the gas phase.
//   const std::vector<Phase> phases_;

//   micm::SystemParameters parameters;

//   parameters.gas_phase_ = micm::Phase();

//   parameters.phases_ = std::vector<micm::Phase> {
//     micm::Phase(),
//     micm::Phase(),
//     micm::Phase(),
  // };

//   return parameters;
// }

// TEST(System, ConstructorWithPhase){
//   micm::SystemParameters parameters;

//   parameters.gas_phase_ = micm::Phase();

//   micm::System system{parameters};
//   system.parameters_.gas_phase_ = micm::Phase();

//   EXPECT_EQ(system.gas_phase_.species_.size(), 0);
// }

// TEST(System, ConstructorWithMultiplePhases){
//   micm::SystemParameters parameters;
//   parameters.phases_ = std::vector<micm::Phase> {
//     micm::Phase(),
//     micm::Phase(),
//     micm::Phase(),
//   };

//   micm::System system{parameters};

//   EXPECT_EQ(system.phases_.size(), 3);
// }

// TEST(System, ConstructorWithAllParameters){
//   micm::System system{FullSetOfParameters()};

//   EXPECT_EQ(system.gas_phase_.species_.size(), 0);
//   EXPECT_EQ(system.phases_.size(), 3);
// }

TEST(System, ConstructorWithParameters){
  std::vector<micm::Species> speciesA = {micm::Species("species1"), micm::Species("species2")};
  std::vector<micm::Species> speciesB = {micm::Species("species3"), micm::Species("species4")};
    
  micm::Phase phase = speciesA;  
  std::vector<micm::Phase> phases = {speciesA, speciesB};
    
  micm::System system = {phase, phases};

  EXPECT_EQ(system.gas_phase_.species_.size(), 2);
  EXPECT_EQ(system.phases_.size(), 2);
}