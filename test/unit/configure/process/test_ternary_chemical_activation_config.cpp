#include <gtest/gtest.h>

#include <micm/configure/solver_config.hpp>

TEST(TernaryChemicalActivationConfig, DetectsInvalidConfig)
{
  micm::SolverConfig solverConfig;

  // Read and parse the configure files
  micm::ConfigParseStatus status = solverConfig.ReadAndParse("./unit_configs/process/ternary_chemical_activation/missing_reactants");
  EXPECT_EQ(micm::ConfigParseStatus::RequiredKeyNotFound, status);
}

