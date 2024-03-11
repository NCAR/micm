#include <gtest/gtest.h>

#include <micm/configure/solver_config.hpp>

TEST(SolverConfig, DetectsInvalidConfigFile)
{
  micm::SolverConfig solverConfig{};
  auto status = solverConfig.ReadAndParse("not_a_config_file");
  EXPECT_EQ(micm::ConfigParseStatus::InvalidCAMPFilePath, status);
}

TEST(SolverConfig, NoConfigFilesFound)
{
  micm::SolverConfig solverConfig{};
  auto status = solverConfig.ReadAndParse("./unit_configs/CAMP/camp_invalid/config.json");
  EXPECT_EQ(micm::ConfigParseStatus::NoConfigFilesFound, status);
}

TEST(SolverConfig, UnknownKey)
{
  micm::SolverConfig solverConfig{};
  auto status = solverConfig.ReadAndParse("./unit_configs/CAMP/camp_unknown_key/config.json");
  EXPECT_EQ(micm::ConfigParseStatus::UnknownKey, status);
  try
  {
    micm::SolverParameters solver_params = solverConfig.GetSolverParams();
  }
  catch (std::runtime_error)
  {
  }
}

TEST(SolverConfig, BadType)
{
  micm::SolverConfig solverConfig{};
  auto status = solverConfig.ReadAndParse("./unit_configs/CAMP/camp_bad_type/config.json");
  EXPECT_EQ(micm::ConfigParseStatus::UnknownKey, status);
  try
  {
    micm::SolverParameters solver_params = solverConfig.GetSolverParams();
  }
  catch (std::runtime_error)
  {
  }
}

TEST(SolverConfig, ReadAndParseCAMPFiles)
{
  micm::SolverConfig solverConfig{};

  // Read and parse the CAMP configure file
  micm::ConfigParseStatus status = solverConfig.ReadAndParse("./unit_configs/CAMP/camp_valid/config.json");
  EXPECT_EQ(micm::ConfigParseStatus::Success, status);
}

TEST(SolverConfig, ReadAndParseCAMPFilesFromDir)
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
  micm::ConfigParseStatus status = solverConfig.ReadAndParse("./unit_configs/chapman/config.json");
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

TEST(SolverConfig, ReadAndParseProcessObjects)
{
  micm::SolverConfig<micm::JsonReaderPolicy> solverConfig;

  // Read and parse the configure files
  micm::ConfigParseStatus status = solverConfig.ReadAndParse("./unit_configs/chapman/config.json");
  EXPECT_EQ(micm::ConfigParseStatus::Success, status);

  // Get solver parameters ('System', the collection of 'Process')
  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  auto& process_vector = solver_params.processes_;

  // Check the number of 'Process' created
  EXPECT_EQ(process_vector.size(), 15);

  // Check the number of 'reactants' and 'products' in each 'Process'
  // Check 'yield' value for the first product and the number of 'spieces in 'phase' in each 'Process'
  int num_reactants_in_each_process[] = { 1, 1, 1, 2, 2, 2, 3, 0, 0, 0, 1, 1, 1, 1, 1 };
  int num_products_in_each_process[] = { 1, 2, 2, 2, 2, 1, 2, 1, 1, 1, 0, 0, 0, 0, 0 };
  double yield_value_of_first_product_in_each_process[] = { 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0,
                                                            1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  int num_phase_in_each_process = 9;

  short idx = 0;
  for (const auto& p : process_vector)
  {
    EXPECT_EQ(p.reactants_.size(), num_reactants_in_each_process[idx]);
    EXPECT_EQ(p.products_.size(), num_products_in_each_process[idx]);
    if (num_products_in_each_process[idx] > 0)
      EXPECT_EQ(p.products_[0].second, yield_value_of_first_product_in_each_process[idx]);
    EXPECT_EQ(p.phase_.species_.size(), num_phase_in_each_process);
    idx++;
  }

  // Check the parameters for 'ArrheniusRateConstant'
  micm::ArrheniusRateConstant* arrhenius_rate_const = nullptr;
  double A_param[] = { 2.15e-11, 3.3e-11, 8.0e-12, 6.0e-34 };
  double B_param[] = { 0.0, 0.0, 0.0, 2.4 };
  double C_param[] = { 110.0, 55.0, -2060.00, 0.0 };
  double D_param[] = { 300.0, 300.0, 300.0, 300.0 };
  double E_param[] = { 0.0, 0.0, 0.0, 0.0 };

  // Convert Arrhenius parameters from expecting molecules cm-3 to moles m-3
  const double conv = 1.0e-6 * 6.02214076e23;
  A_param[0] *= conv;         // 2 reactants
  A_param[1] *= conv;         // 2 reactants
  A_param[2] *= conv;         // 2 reactants
  A_param[3] *= conv * conv;  // 3 reactants

  for (short i = 3; i < 7; i++)
  {
    arrhenius_rate_const = dynamic_cast<micm::ArrheniusRateConstant*>(process_vector[i].rate_constant_.get());

    EXPECT_TRUE(arrhenius_rate_const != nullptr);
    EXPECT_EQ(arrhenius_rate_const->parameters_.A_, A_param[i - 3]);
    EXPECT_EQ(arrhenius_rate_const->parameters_.B_, B_param[i - 3]);
    EXPECT_EQ(arrhenius_rate_const->parameters_.C_, C_param[i - 3]);
    EXPECT_EQ(arrhenius_rate_const->parameters_.D_, D_param[i - 3]);
    EXPECT_EQ(arrhenius_rate_const->parameters_.E_, E_param[i - 3]);
  }

  // Check the number of custom parameters of 'rate constant' in each 'Process'
  std::vector<std::vector<std::string>> custom_rate_labels = { { "PHOTO.O2_1" },
                                                               { "PHOTO.O3_1" },
                                                               { "PHOTO.O3_2" },
                                                               {},
                                                               {},
                                                               {},
                                                               {},
                                                               { "EMIS.O1D" },
                                                               { "EMIS.O" },
                                                               { "EMIS.O3" },
                                                               { "LOSS.N2" },
                                                               { "LOSS.O2" },
                                                               { "LOSS.CO2" },
                                                               { "LOSS.Ar" },
                                                               { "LOSS.H2O" } };

  // check photlysis, emissions, and loss reaction labels
  idx = 0;
  std::vector<micm::Process>::iterator it;
  for (it = solver_params.processes_.begin(); it != solver_params.processes_.end(); it++, idx++)
  {
    EXPECT_EQ(it->rate_constant_->SizeCustomParameters(), custom_rate_labels[idx].size());
    for (std::size_t i = 0; i < custom_rate_labels[idx].size(); ++i)
      EXPECT_EQ(it->rate_constant_->CustomParameters()[i], custom_rate_labels[idx][i]);
  }
}
