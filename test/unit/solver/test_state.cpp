#include <gtest/gtest.h>

#include <micm/process/photolysis_rate_constant.hpp>
#include <micm/solver/state.hpp>

TEST(State, Constructor)
{
  micm::State<micm::Matrix> state{ micm::StateParameters{ .state_variable_names_{ "foo", "bar", "baz", "quz" },
                                                          .number_of_grid_cells_ = 3,
                                                          .number_of_custom_parameters_ = 5,
                                                          .number_of_rate_constants_ = 10 } };

  EXPECT_EQ(state.conditions_.size(), 3);
  EXPECT_EQ(state.variable_map_["foo"], 0);
  EXPECT_EQ(state.variable_map_["bar"], 1);
  EXPECT_EQ(state.variable_map_["baz"], 2);
  EXPECT_EQ(state.variable_map_["quz"], 3);
  EXPECT_EQ(state.variables_.size(), 3);
  EXPECT_EQ(state.variables_[0].size(), 4);
  EXPECT_EQ(state.custom_rate_parameters_.size(), 3);
  EXPECT_EQ(state.custom_rate_parameters_[0].size(), 5);
  EXPECT_EQ(state.rate_constants_.size(), 3);
  EXPECT_EQ(state.rate_constants_[0].size(), 10);
}

TEST(State, SettingConcentrationsWithInvalidArguementsThrowsException)
{
  micm::State<micm::Matrix> state{ micm::StateParameters{ .state_variable_names_{ "foo", "bar", "baz", "quz" },
                                                          .number_of_grid_cells_ = 3,
                                                          .number_of_custom_parameters_ = 5,
                                                          .number_of_rate_constants_ = 10 } };

  std::unordered_map<std::string, std::vector<double>> concentrations = {
    { "FUU", { 0.1 } }, { "bar", { 0.2 } }, { "baz", { 0.3 } }, { "quz", { 0.4 } }
  };

  // Build system
  micm::Phase gas_phase(
      std::vector<micm::Species>{ micm::Species("foo"), micm::Species("bar"), micm::Species("baz"), micm::Species("quz") });
  micm::System system{ micm::SystemParameters{ gas_phase } };

  EXPECT_ANY_THROW(state.SetConcentrations(system, concentrations));
}

TEST(State, SetConcentrations)
{
  uint32_t num_grid_cells = 3;
  uint32_t num_species = 4;

  micm::State<micm::Matrix> state{ micm::StateParameters{ .state_variable_names_{ "foo", "bar", "baz", "quz" },
                                                          .number_of_grid_cells_ = num_grid_cells,
                                                          .number_of_custom_parameters_ = 5,
                                                          .number_of_rate_constants_ = 10 } };

  // Build system
  micm::Phase gas_phase(
      std::vector<micm::Species>{ micm::Species("foo"), micm::Species("bar"), micm::Species("baz"), micm::Species("quz") });
  micm::System system{ micm::SystemParameters{ gas_phase } };

  std::unordered_map<std::string, std::vector<double>> concentrations = { { "bar", { 0.2, 0.22, 0.222 } },
                                                                          { "baz", { 0.3, 0.33, 0.333 } },
                                                                          { "foo", { 0.9, 0.99, 0.999 } },
                                                                          { "quz", { 0.4, 0.44, 0.444 } } };

  state.SetConcentrations(system, concentrations);

  // Compare concentration values
  std::vector<double> concentrations_in_order{
    0.9, 0.2, 0.3, 0.4, 0.99, 0.22, 0.33, 0.44, 0.999, 0.222, 0.333, 0.444,
  };

  short idx = 0;
  for (short i = 0; i < num_grid_cells; i++)
  {
    for (short j = 0; j < num_species; j++, idx++)
    {
      EXPECT_EQ(state.variables_[i][j], concentrations_in_order[idx]);
    }
  }
}

TEST(State, SettingPhotolysisRateWithInvalidArguementsThrowsException)
{
  micm::State<micm::Matrix> state{ micm::StateParameters{ .state_variable_names_{ "foo", "bar", "baz", "quz" },
                                                          .number_of_grid_cells_ = 3,
                                                          .number_of_custom_parameters_ = 5,
                                                          .number_of_rate_constants_ = 10 } };

  // Build an array of photolysis rate constant
  std::vector<micm::PhotolysisRateConstant> photolysis_rate_arr;
  photolysis_rate_arr.reserve(5);
  photolysis_rate_arr.emplace_back("O1");
  photolysis_rate_arr.emplace_back("O2");
  photolysis_rate_arr.emplace_back("O3");
  photolysis_rate_arr.emplace_back("AAA");
  photolysis_rate_arr.emplace_back("BBB");

  // user input for photolysis rate constant (unordered)
  std::unordered_map<std::string, std::vector<double>> photo_rates = {
    { "O3", { 0.3 } }, { "O1", { 0.1 } }, { "O2", { 0.5 } }, { "CCC", { 0.7 } }, { "AAA", { 0.5 } }
  };

  std::vector<double> photo_rates_in_order{ 0.1, 0.5, 0.3, 0.5, 0.7 };

  // Compare photolysis rate constants
  EXPECT_ANY_THROW(state.SetPhotolysisRate(photolysis_rate_arr, photo_rates));
}

TEST(State, SetPhotolysisRate)
{
  uint32_t num_grid_cells = 3;
  uint32_t num_custom_params = 5;

  micm::State<micm::Matrix> state{ micm::StateParameters{ .state_variable_names_{ "foo", "bar", "baz", "quz" },
                                                          .number_of_grid_cells_ = num_grid_cells,
                                                          .number_of_custom_parameters_ = num_custom_params,
                                                          .number_of_rate_constants_ = 10 } };

  // Build an array of photolysis rate constant
  std::vector<micm::PhotolysisRateConstant> photolysis_rate_arr;
  photolysis_rate_arr.reserve(num_custom_params);
  photolysis_rate_arr.emplace_back("O1");
  photolysis_rate_arr.emplace_back("O2");
  photolysis_rate_arr.emplace_back("O3");
  photolysis_rate_arr.emplace_back("AAA");
  photolysis_rate_arr.emplace_back("BBB");

  // user input for photolysis rate constant (unordered)
  std::unordered_map<std::string, std::vector<double>> photo_rates = { { "O3", { 0.3, 0.33, 0.333 } },
                                                                       { "O1", { 0.1, 0.11, 0.111 } },
                                                                       { "O2", { 0.5, 0.55, 0.555 } },
                                                                       { "BBB", { 0.7, 0.77, 0.777 } },
                                                                       { "AAA", { 0.5, 0.55, 0.555 } } };

  std::vector<double> photo_rates_in_order{ 0.1,  0.5,  0.3,   0.5,   0.7,   0.11,  0.55, 0.33,
                                            0.55, 0.77, 0.111, 0.555, 0.333, 0.555, 0.777 };

  // Compare photolysis rate constants
  state.SetPhotolysisRate(photolysis_rate_arr, photo_rates);

  short idx = 0;
  for (short i = 0; i < num_grid_cells; i++)
  {
    for (short j = 0; j < num_custom_params; j++, idx++)
    {
      EXPECT_EQ(state.custom_rate_parameters_[i][j], photo_rates_in_order[idx]);
    }
  }
}
