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

  std::unordered_map<std::string, double> concentrations = {
    { "FUU", 0.1 }, { "bar", 0.2 }, { "baz", 0.3 }, { "quz", 0.4 }
  };

  // Build system
  micm::Phase gas_phase(
      std::vector<micm::Species>{ micm::Species("foo"), micm::Species("bar"), micm::Species("baz"), micm::Species("quz") });
  micm::System system{ micm::SystemParameters{ gas_phase } };

  EXPECT_ANY_THROW(state.SetConcentrations(system, concentrations));
}

TEST(State, SetConcentrations)
{
  micm::State<micm::Matrix> state{ micm::StateParameters{ .state_variable_names_{ "foo", "bar", "baz", "quz" },
                                            .number_of_grid_cells_ = 3,
                                            .number_of_custom_parameters_ = 5,
                                            .number_of_rate_constants_ = 10 } };

  // Build system
  micm::Phase gas_phase(
      std::vector<micm::Species>{ micm::Species("foo"), micm::Species("bar"), micm::Species("baz"), micm::Species("quz") });
  micm::System system{ micm::SystemParameters{ gas_phase } };

  std::unordered_map<std::string, double> concentrations = {
    { "bar", 0.2 }, { "baz", 0.3 }, { "foo", 0.99 }, { "quz", 0.4 }
  };

  state.SetConcentrations(system, concentrations);

  // Compare concentration values
  std::vector<double> concentrations_in_order{ 0.99, 0.2, 0.3, 0.4 };

  short idx = 0;
  for (auto& val : state.variables_[0])
  {
    EXPECT_EQ(val, concentrations_in_order[idx]);
    idx++;
  }
}

TEST(State, SetPhotolysisRate)
{
  micm::State<micm::Matrix> state{ micm::StateParameters{ .state_variable_names_{ "foo", "bar", "baz", "quz" },
                                            .number_of_grid_cells_ = 3,
                                            .number_of_custom_parameters_ = 5,
                                            .number_of_rate_constants_ = 10 } };

  // Build an array of photolysis rate constant
  std::vector<micm::PhotolysisRateConstant> photolysis_rate_arr;
  photolysis_rate_arr.reserve(3);
  photolysis_rate_arr.emplace_back("O1");
  photolysis_rate_arr.emplace_back("O2");
  photolysis_rate_arr.emplace_back("O3");

  // user input for photolysis rate constant (unordered)
  std::unordered_map<std::string, double> photo_rates = { { "O3", 0.3 }, { "O1", 0.1 }, { "O2", 0.5 } };

  std::vector<double> photo_rates_in_order{ 0.1, 0.5, 0.3 };

  // Compare photolysis rate constants
  state.SetPhotolysisRate(photolysis_rate_arr, photo_rates);

  short idx = 0;
  for (auto& val : state.custom_rate_parameters_[0])
  {
    EXPECT_EQ(val, photo_rates_in_order[idx]);
    idx++;
  }
}