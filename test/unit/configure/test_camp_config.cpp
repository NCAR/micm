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

TEST(SolverConfig, ReadAndParseSystemObject)
{
  micm::SolverConfig solverConfig;

  // Read and parse the configure files
  micm::ConfigParseStatus status = solverConfig.ReadAndParse("./unit_configs/chapman");
  EXPECT_EQ(micm::ConfigParseStatus::Success, status);

  // Get solver parameters ('System', the collection of 'Process')
  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  // Check 'gas_phase' in 'System'
  EXPECT_EQ(solver_params.system_.gas_phase_.species_.size(), 9);

  // Check 'phases' in 'System'
  EXPECT_EQ(solver_params.system_.phases_.size(), 0);

  // Check 'name' and 'properties' in 'Species'
  std::vector<std::pair<std::string, short>> species_name_and_num_properties = {
    std::make_pair("M", 0),   std::make_pair("Ar", 1), std::make_pair("CO2", 1),
    std::make_pair("H2O", 1), std::make_pair("N2", 1), std::make_pair("O1D", 1),
    std::make_pair("O", 1),   std::make_pair("O2", 1), std::make_pair("O3", 1)
  };

  short idx = 0;
  for (const auto& s : solver_params.system_.gas_phase_.species_)
  {
    EXPECT_EQ(s.name_, species_name_and_num_properties[idx].first);
    EXPECT_EQ(s.properties_.size(), species_name_and_num_properties[idx].second);
    idx++;
  }
}
