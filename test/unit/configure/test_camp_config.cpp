#include <gtest/gtest.h>

#include <micm/configure/camp_config.hpp>

/*
TEST(SolverConfig, DetectsInvalidConfigFile)
{
  micm::SolverConfig solverConfig{};
  auto status = solverConfig.ReadAndParse("not_a_config_file_directory");
  EXPECT_EQ(micm::ConfigParseStatus::InvalidCAMPFilePath, status);
}
*/

TEST(SolverConfig, ReadAndParseCAMPFiles)
{
  micm::SolverConfig solverConfig{};

  // Read and parse the CAMP configure file
  micm::ConfigParseStatus status = solverConfig.ReadAndParse("./unit_configs/CAMP/camp_valid");
  EXPECT_EQ(micm::ConfigParseStatus::Success, status);
}

