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
  auto* solver_params_ptr = std::get_if<micm::SolverParameters>(&configs);

  EXPECT_TRUE(solver_params_ptr != nullptr);

  micm::SolverParameters& solver_params = *solver_params_ptr;

  EXPECT_EQ(solver_params.system_.gas_phase_.species_.size(), 9);
  EXPECT_EQ(solver_params.system_.phases_.size(), 0);

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
  auto* solver_params_ptr = std::get_if<micm::SolverParameters>(&configs);

  EXPECT_TRUE(solver_params_ptr != nullptr);

  micm::SolverParameters& solver_params = *solver_params_ptr;

  int num_processes = 7;
  int num_reactants_in_each_process[] = { 1, 1, 1, 2, 2, 2, 3 };
  int num_products_in_each_process[] = { 1, 2, 2, 2, 2, 1, 2 };
  double yield_value_of_first_product_in_each_process[] = { 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0 };
  std::size_t size_custom_parameters_of_rate_constant_in_each_process[] = { 1, 1, 1, 0, 0, 0, 0 };
  int num_phase_in_each_process = 9;

  EXPECT_EQ(solver_params.processes_.size(), num_processes);

  short idx = 0;
  for (const auto& p : solver_params.processes_)
  {
    EXPECT_EQ(p.reactants_.size(), num_reactants_in_each_process[idx]);
    EXPECT_EQ(p.products_.size(), num_products_in_each_process[idx]);
    EXPECT_EQ(p.products_[0].second, yield_value_of_first_product_in_each_process[idx]);
    EXPECT_EQ(p.phase_.species_.size(), num_phase_in_each_process);
    idx++;
  }

  idx = 0;
  std::vector<micm::Process>::iterator it;
  for (it = solver_params.processes_.begin(); it != solver_params.processes_.end(); it++, idx++)
  {
    EXPECT_EQ(it->rate_constant_->SizeCustomParameters(), size_custom_parameters_of_rate_constant_in_each_process[idx]);
  }
}
#endif