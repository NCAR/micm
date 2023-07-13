#include <gtest/gtest.h>

#ifdef USE_JSON
#  include <micm/configure/solver_config.hpp>

TEST(SolverConfig, DetectsInvalidConfigFileAndThrow)
{
  micm::SolverConfig<micm::JsonReaderPolicy, micm::ThrowPolicy> solverConfig{};
  EXPECT_ANY_THROW(solverConfig.Configure("not_a_config_file.json"));
}

TEST(SolverConfig, DetectsInvalidConfigFileAndNoThrowDoesntThrow)
{
  micm::SolverConfig<micm::JsonReaderPolicy, micm::NoThrowPolicy> solverConfig{};
  EXPECT_NO_THROW(solverConfig.Configure("not_a_config_file.json"));

  std::variant<micm::SolverParameters, micm::ConfigErrorCode> configs = solverConfig.Configure("not_a_config_file.json");
  EXPECT_EQ(std::get<micm::ConfigErrorCode>(configs), micm::ConfigErrorCode::FileNotFound);
}

TEST(SolverConfig, ReadAndParseSystem)
{
  micm::SolverConfig<micm::JsonReaderPolicy, micm::ThrowPolicy> solverConfig{};
  std::variant<micm::SolverParameters, micm::ConfigErrorCode> configs =
      solverConfig.Configure("./unit_configs/chapman/config.json");

  // Check if parsing is successful and returns 'Solverparameters'
  auto* solver_params_ptr = std::get_if<micm::SolverParameters>(&configs);
  EXPECT_TRUE(solver_params_ptr != nullptr);

  micm::SolverParameters& solver_params = *solver_params_ptr;

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

TEST(SolverConfig, ReadAndParseProcess)
{
  micm::SolverConfig<micm::JsonReaderPolicy, micm::ThrowPolicy> solverConfig{};
  std::variant<micm::SolverParameters, micm::ConfigErrorCode> configs =
      solverConfig.Configure("./unit_configs/chapman/config.json");

  // Check if parsing is successful and returns 'Solverparameters'
  auto* solver_params_ptr = std::get_if<micm::SolverParameters>(&configs);
  EXPECT_TRUE(solver_params_ptr != nullptr);

  micm::SolverParameters& solver_params = *solver_params_ptr;
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
#endif