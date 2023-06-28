#include <gtest/gtest.h>

#include <micm/solver/state.hpp>

TEST(State, Constructor)
{
  micm::State<micm::Matrix> state {micm::StateParameters{
    .state_variable_names_{ "foo", "bar", "baz", "quz" },
    .number_of_grid_cells_ = 3,
    .number_of_custom_parameters_ = 5,
    .number_of_rate_constants_ = 10
  }};

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
  micm::State<micm::Matrix> state {micm::StateParameters{
    .state_variable_names_{ "foo", "bar", "baz", "quz" },
    .number_of_grid_cells_ = 3,
    .number_of_custom_parameters_ = 5,
    .number_of_rate_constants_ = 10
  }};

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
  micm::State state{ micm::StateParameters{ .state_variable_names_{ "foo", "bar", "baz", "quz" },
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
