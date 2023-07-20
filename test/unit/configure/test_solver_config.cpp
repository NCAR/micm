#include <gtest/gtest.h>

#ifdef USE_JSON
#  include <micm/configure/solver_config.hpp>

TEST(SolverConfig, DetectsInvalidConfigFileAndThrow)
{
  micm::SolverConfig<micm::JsonReaderPolicy, micm::ThrowPolicy> solverConfig{};
  EXPECT_ANY_THROW(solverConfig.ReadAndParse("not_a_config_file_directory"));
}

TEST(SolverConfig, DetectsInvalidConfigFileAndNoThrowDoesntThrow)
{
  micm::SolverConfig<micm::JsonReaderPolicy, micm::NoThrowPolicy> solverConfig{};
  EXPECT_NO_THROW(solverConfig.ReadAndParse("not_a_config_file_directory"));
  EXPECT_FALSE(solverConfig.ReadAndParse("not_a_config_file_directory"));
}

TEST(SolverConfig, ReadAndParseSystemObject)
{
  micm::SolverConfig<micm::JsonReaderPolicy, micm::ThrowPolicy> solverConfig;  // Set to throw-exception policy

  // Read and parse the configure files
  bool is_parse_success = solverConfig.ReadAndParse("./unit_configs/chapman");
  EXPECT_TRUE(is_parse_success);

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
  micm::SolverConfig<micm::JsonReaderPolicy, micm::ThrowPolicy> solverConfig;  // Set to throw-exception policy

  // Read and parse the configure files
  bool is_parse_success = solverConfig.ReadAndParse("./unit_configs/chapman");
  EXPECT_TRUE(is_parse_success);

  // Get solver parameters ('System', the collection of 'Process')
  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  auto& process_vector = solver_params.processes_;

  // Check the number of 'Process' created
  EXPECT_EQ(process_vector.size(), 7);

  // Check the number of 'reactants' and 'products' in each 'Process'
  // Check 'yield' value for the first product and the number of 'spieces in 'phase' in each 'Process'
  int num_reactants_in_each_process[] = { 1, 1, 1, 2, 2, 2, 3 };
  int num_products_in_each_process[] = { 1, 2, 2, 2, 2, 1, 2 };
  double yield_value_of_first_product_in_each_process[] = { 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0 };
  int num_phase_in_each_process = 9;

  short idx = 0;
  for (const auto& p : process_vector)
  {
    EXPECT_EQ(p.reactants_.size(), num_reactants_in_each_process[idx]);
    EXPECT_EQ(p.products_.size(), num_products_in_each_process[idx]);
    EXPECT_EQ(p.products_[0].second, yield_value_of_first_product_in_each_process[idx]);
    EXPECT_EQ(p.phase_.species_.size(), num_phase_in_each_process);
    idx++;
  }

  // Check the name for 'PhotolysisRateConstant'
  micm::PhotolysisRateConstant* photolysis_rate_const = nullptr;
  std::string photolysis_name[] = { "O2_1", "O3_1", "O3_2" };

  for (short i = 0; i < 3; i++)
  {
    photolysis_rate_const = dynamic_cast<micm::PhotolysisRateConstant*>(process_vector[i].rate_constant_.get());

    EXPECT_TRUE(photolysis_rate_const != nullptr);
    EXPECT_EQ(photolysis_rate_const->name_, photolysis_name[i]);
  }

  // Check the parameters for 'ArrheniusRateConstant'
  micm::ArrheniusRateConstant* arrhenius_rate_const = nullptr;
  double A_param[] = { 2.15e-11, 3.3e-11, 8.0e-12, 6.0e-34 };
  double B_param[] = { 0.0, 0.0, 0.0, 2.4 };
  double C_param[] = { 110.0, 55.0, -2060.00, 0.0 };
  double D_param[] = { 300.0, 300.0, 300.0, 300.0 };
  double E_param[] = { 0.0, 0.0, 0.0, 0.0 };

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
  std::size_t size_custom_parameters_of_rate_constant_in_each_process[] = { 1, 1, 1, 0, 0, 0, 0 };

  idx = 0;
  std::vector<micm::Process>::iterator it;
  for (it = solver_params.processes_.begin(); it != solver_params.processes_.end(); it++, idx++)
  {
    EXPECT_EQ(it->rate_constant_->SizeCustomParameters(), size_custom_parameters_of_rate_constant_in_each_process[idx]);
  }
}

TEST(SolverConfig, GettingSolverParamsThrowsExceptionWithFailedParsing)
{
  micm::SolverConfig<micm::JsonReaderPolicy, micm::NoThrowPolicy> solverConfig;
  bool is_parse_success = solverConfig.ReadAndParse("not_a_config_file_directory");
  EXPECT_FALSE(is_parse_success);
  EXPECT_ANY_THROW(solverConfig.GetSolverParams());
}

TEST(SolverConfig, GettingPhotolysisRateConstantThrowsExceptionWithFailedParsing)
{
  micm::SolverConfig<micm::JsonReaderPolicy, micm::NoThrowPolicy> solverConfig;
  bool is_parse_success = solverConfig.ReadAndParse("not_a_config_file_directory");
  EXPECT_FALSE(is_parse_success);
  EXPECT_ANY_THROW(solverConfig.GetPhotolysisRateConstants());
}

TEST(SolverConfig, GetPhotolysisRateConstants)
{
  // Read and parse the configure files
  micm::SolverConfig<micm::JsonReaderPolicy, micm::NoThrowPolicy> solverConfig;
  bool is_parse_success = solverConfig.ReadAndParse("./unit_configs/chapman");
  EXPECT_TRUE(is_parse_success);

  std::vector<micm::PhotolysisRateConstant>& photolysis_rate_arr_ = solverConfig.GetPhotolysisRateConstants();

  // Check the name of photolsis rate constants
  std::array<std::string, 3> photo_names{ "O2_1", "O3_1", "O3_2" };

  short idx = 0;
  for (auto& rate : photolysis_rate_arr_)
  {
    EXPECT_EQ(rate.name_, photo_names[idx]);
    idx++;
  }
}

//
// Tests for MZ326 configure files
//
TEST(SolverConfig, ReadAndParseSystemObjectfromMZ326)
{
  micm::SolverConfig<micm::JsonReaderPolicy, micm::ThrowPolicy> solverConfig;  // Set to throw-exception policy

  // Read and parse the configure files
  bool is_parse_success = solverConfig.ReadAndParse("./unit_configs/MZ326");
  EXPECT_TRUE(is_parse_success);

  // Get solver parameters ('System', the collection of 'Process')
  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  // Check 'gas_phase' in 'System'
  EXPECT_EQ(solver_params.system_.gas_phase_.species_.size(), 6);

  // Check 'phases' in 'System'
  EXPECT_EQ(solver_params.system_.phases_.size(), 0);

  // Check 'name' and molecular_weight from 'properties' in 'Species'
  std::vector<std::pair<std::string, double>> species_name_and_molecular_weight = {
    std::make_pair("ALKNIT", 0.133141), std::make_pair("BZOOH", 0.124135), std::make_pair("C6H5OOH", 0.110109),
    std::make_pair("COF2", 0.0),        std::make_pair("O2", 0.0),         std::make_pair("FUR2O2", 0.0)
  };

  short idx = 0;
  for (const auto& s : solver_params.system_.gas_phase_.species_)
  {
    EXPECT_EQ(s.name_, species_name_and_molecular_weight[idx].first);
    EXPECT_EQ(s.properties_[0].value_, species_name_and_molecular_weight[idx].second);
    idx++;
  }
}

TEST(SolverConfig, ReadAndParseProcessObjectsfromMZ326)
{
  micm::SolverConfig<micm::JsonReaderPolicy, micm::ThrowPolicy> solverConfig;  // Set to throw-exception policy

  // Read and parse the configure files
  bool is_parse_success = solverConfig.ReadAndParse("./unit_configs/MZ326");
  EXPECT_TRUE(is_parse_success);

  // Get solver parameters ('System', the collection of 'Process')
  micm::SolverParameters solver_params = solverConfig.GetSolverParams();

  auto& process_vector = solver_params.processes_;

  // Check the number of 'Process' created
  EXPECT_EQ(process_vector.size(), 5);

  // Check the number of 'reactants' and 'products' in each 'Process'
  // Check 'yield' value for the first product and the number of 'spieces in 'phase' in each 'Process'
  int num_reactants_in_each_process[] = { 2, 2, 2, 1, 1 };
  int num_products_in_each_process[] = { 5, 5, 2, 3, 2 };
  int num_phase_in_each_process = 6;

  short idx = 0;
  for (const auto& p : process_vector)
  {
    EXPECT_EQ(p.reactants_.size(), num_reactants_in_each_process[idx]);
    EXPECT_EQ(p.products_.size(), num_products_in_each_process[idx]);
    EXPECT_EQ(p.phase_.species_.size(), num_phase_in_each_process);
    idx++;
  }

  // Check the parameters for 'TroeRateConstant'
  micm::TroeRateConstant* troe_rate_const = nullptr;
  double k0_A_param[] = { 5.5e-30, 1.07767 };
  double k0_B_param[] = { 0.0, -5.6 };
  double k0_C_param[] = { 0.0, -14000 };
  double kinf_A_param[] = { 8.3e-13, 1.03323e+17 };
  double kinf_B_param[] = { 0.0, 0.0 };
  double kinf_C_param[] = { 0.0, -14000 };
  double Fc_param[] = { 0.6, 0.6 };
  double N_param[] = { -2, 1.5 };

  idx = 0;
  for (short i : { 0, 4 })
  {
    troe_rate_const = dynamic_cast<micm::TroeRateConstant*>(process_vector[i].rate_constant_.get());

    EXPECT_TRUE(troe_rate_const != nullptr);
    EXPECT_EQ(troe_rate_const->parameters_.k0_A_, k0_A_param[idx]);
    EXPECT_EQ(troe_rate_const->parameters_.k0_B_, k0_B_param[idx]);
    EXPECT_EQ(troe_rate_const->parameters_.k0_C_, k0_C_param[idx]);
    EXPECT_EQ(troe_rate_const->parameters_.kinf_A_, kinf_A_param[idx]);
    EXPECT_EQ(troe_rate_const->parameters_.kinf_B_, kinf_B_param[idx]);
    EXPECT_EQ(troe_rate_const->parameters_.kinf_C_, kinf_C_param[idx]);
    EXPECT_EQ(troe_rate_const->parameters_.Fc_, Fc_param[idx]);
    EXPECT_EQ(troe_rate_const->parameters_.N_, N_param[idx]);

    idx++;
  }

  // Check the name for 'PhotolysisRateConstant'
  micm::PhotolysisRateConstant* photolysis_rate_const = nullptr;
  std::string photolysis_name = "jterpnit";

  for (short i : { 3 })
  {
    photolysis_rate_const = dynamic_cast<micm::PhotolysisRateConstant*>(process_vector[i].rate_constant_.get());

    EXPECT_TRUE(photolysis_rate_const != nullptr);
    EXPECT_EQ(photolysis_rate_const->name_, photolysis_name);
  }

  // Check the parameters for 'ArrheniusRateConstant'
  micm::ArrheniusRateConstant* arrhenius_rate_const = nullptr;
  double A_param[] = { 2.0e-12, 3.8e-12 };
  double B_param[] = { 0.0, 0.0 };
  double C_param[] = { -1 * -6.90325e-21 / 1.3806505e-23, -1 * -2.7613e-21 / 1.3806505e-23 };
  double D_param[] = { 300.0, 300.0 };
  double E_param[] = { 0.0, 0.0 };

  idx = 0;
  for (short i : { 1, 2 })
  {
    arrhenius_rate_const = dynamic_cast<micm::ArrheniusRateConstant*>(process_vector[i].rate_constant_.get());

    EXPECT_TRUE(arrhenius_rate_const != nullptr);
    EXPECT_DOUBLE_EQ(arrhenius_rate_const->parameters_.A_, A_param[idx]);
    EXPECT_DOUBLE_EQ(arrhenius_rate_const->parameters_.B_, B_param[idx]);
    EXPECT_DOUBLE_EQ(arrhenius_rate_const->parameters_.C_, C_param[idx]);
    EXPECT_DOUBLE_EQ(arrhenius_rate_const->parameters_.D_, D_param[idx]);
    EXPECT_DOUBLE_EQ(arrhenius_rate_const->parameters_.E_, E_param[idx]);

    idx++;
  }
}
#endif
