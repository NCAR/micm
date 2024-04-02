#include <gtest/gtest.h>

#include <micm/solver/state.hpp>

TEST(State, Constructor)
{
  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = 3,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "quux", "corge" },
  } };

  EXPECT_EQ(state.conditions_.size(), 3);
  EXPECT_EQ(state.variable_map_["foo"], 0);
  EXPECT_EQ(state.variable_map_["bar"], 1);
  EXPECT_EQ(state.variable_map_["baz"], 2);
  EXPECT_EQ(state.variable_map_["quz"], 3);
  EXPECT_EQ(state.variables_.size(), 3);
  EXPECT_EQ(state.variables_[0].size(), 4);
  EXPECT_EQ(state.custom_rate_parameter_map_["quux"], 0);
  EXPECT_EQ(state.custom_rate_parameter_map_["corge"], 1);
  EXPECT_EQ(state.custom_rate_parameters_.size(), 3);
  EXPECT_EQ(state.custom_rate_parameters_[0].size(), 2);
  EXPECT_EQ(state.rate_constants_.size(), 3);
  EXPECT_EQ(state.rate_constants_[0].size(), 10);
}

TEST(State, SettingSingleConcentrationWithInvalidArgumentsThowsException)
{
  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = 3,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "quux", "corge" },
  } };
  EXPECT_ANY_THROW(state.SetConcentration(micm::Species{ "foo" }, 1.0));
  EXPECT_ANY_THROW(state.SetConcentration(micm::Species{ "foo" }, std::vector<double>{ 1.0, 2.0 }));
  EXPECT_ANY_THROW(state.SetConcentration(micm::Species{ "not foo" }, 1.0));
  EXPECT_ANY_THROW(state.SetConcentration(micm::Species{ "not foo" }, std::vector<double>{ 1.0, 2.0, 3.0, 4.0 }));
}

TEST(State, SetSingleConcentration)
{
  {
    micm::State state{ micm::StateParameters{
        .number_of_grid_cells_ = 3,
        .number_of_rate_constants_ = 10,
        .variable_names_{ "foo", "bar", "baz", "quz" },
        .custom_rate_parameter_labels_{ "quux", "corge" },
    } };
    std::vector<double> concentrations{ 12.0, 42.0, 35.2 };
    state.SetConcentration(micm::Species{ "bar" }, concentrations);
    for (std::size_t i = 0; i < concentrations.size(); ++i)
      EXPECT_EQ(state.variables_[i][state.variable_map_["bar"]], concentrations[i]);
  }
  {
    micm::State state{ micm::StateParameters{
        .number_of_grid_cells_ = 1,
        .number_of_rate_constants_ = 10,
        .variable_names_{ "foo", "bar", "baz", "quz" },
        .custom_rate_parameter_labels_{ "quux", "corge" },
    } };
    state.SetConcentration(micm::Species{ "bar" }, 324.2);
    EXPECT_EQ(state.variables_[0][state.variable_map_["bar"]], 324.2);
  }
}

TEST(State, SettingConcentrationsWithInvalidArguementsThrowsException)
{
  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = 3,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "quux", "corge" },
  } };

  std::unordered_map<std::string, std::vector<double>> concentrations = {
    { "FUU", { 0.1 } }, { "bar", { 0.2 } }, { "baz", { 0.3 } }, { "quz", { 0.4 } }
  };

  EXPECT_ANY_THROW(state.SetConcentrations(concentrations));
}

TEST(State, SetConcentrations)
{
  uint32_t num_grid_cells = 3;
  uint32_t num_species = 4;

  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = num_grid_cells,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "quux", "corge" },
  } };

  std::unordered_map<std::string, std::vector<double>> concentrations = { { "bar", { 0.2, 0.22, 0.222 } },
                                                                          { "baz", { 0.3, 0.33, 0.333 } },
                                                                          { "foo", { 0.9, 0.99, 0.999 } },
                                                                          { "quz", { 0.4, 0.44, 0.444 } } };

  state.SetConcentrations(concentrations);

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

TEST(State, SettingCustomRateParameterWithInvalidVectorSizeThrowsException)
{
  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = 3,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
  } };

  // user input for custom rate parameters (unordered)
  std::unordered_map<std::string, std::vector<double>> custom_params = {
    { "O3", { 0.3 } }, { "O1", { 0.1 } }, { "O2", { 0.5 } }, { "BBB", { 0.7 } }, { "AAA", { 0.5 } }
  };

  EXPECT_ANY_THROW(state.SetCustomRateParameters(custom_params));
}

TEST(State, SettingCustomRateParameterWithInvalidLabelThrowsException)
{
  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = 1,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
  } };

  // user input for custom rate parameters (unordered)
  std::unordered_map<std::string, std::vector<double>> custom_params = {
    { "O3", { 0.3 } }, { "O1", { 0.1 } }, { "O2", { 0.5 } }, { "CCC", { 0.7 } }, { "AAA", { 0.5 } }
  };

  EXPECT_ANY_THROW(state.SetCustomRateParameters(custom_params));
}

TEST(State, SetCustomRateParameter)
{
  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = 1,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
  } };

  state.SetCustomRateParameter("O2", 42.3);
  EXPECT_EQ(state.custom_rate_parameters_[0][1], 42.3);
}

TEST(State, SetCustomRateParameters)
{
  uint32_t num_grid_cells = 3;

  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = num_grid_cells,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
  } };

  // user input for custom rate parameters (unordered)
  std::unordered_map<std::string, std::vector<double>> custom_params = { { "O3", { 0.3, 0.33, 0.333 } },
                                                                         { "O1", { 0.1, 0.11, 0.111 } },
                                                                         { "O2", { 0.5, 0.55, 0.555 } },
                                                                         { "BBB", { 0.7, 0.77, 0.777 } },
                                                                         { "AAA", { 0.5, 0.55, 0.555 } } };

  std::vector<double> custom_params_in_order{ 0.1,  0.5,  0.3,   0.5,   0.7,   0.11,  0.55, 0.33,
                                              0.55, 0.77, 0.111, 0.555, 0.333, 0.555, 0.777 };

  state.SetCustomRateParameters(custom_params);

  short idx = 0;
  for (short i = 0; i < num_grid_cells; i++)
  {
    for (short j = 0; j < 5; j++, idx++)
    {
      EXPECT_EQ(state.custom_rate_parameters_[i][j], custom_params_in_order[idx]);
    }
  }
}

TEST(State, UnsafelySetCustomRateParameterOneCell)
{
  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = 1,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
  } };

  std::vector<std::vector<double>> parameters = { { 0.1, 0.2, 0.3, 0.4, 0.5 } };

  state.UnsafelySetCustomRateParameters(parameters);
  EXPECT_EQ(state.custom_rate_parameters_[0][0], 0.1);
  EXPECT_EQ(state.custom_rate_parameters_[0][1], 0.2);
  EXPECT_EQ(state.custom_rate_parameters_[0][2], 0.3);
  EXPECT_EQ(state.custom_rate_parameters_[0][3], 0.4);
  EXPECT_EQ(state.custom_rate_parameters_[0][4], 0.5);
}

TEST(State, UnsafelySetCustomRateParameterMultiCell)
{
  uint32_t num_grid_cells = 3;

  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = num_grid_cells,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
  } };

  std::vector<std::vector<double>> parameters = { { 0.1, 0.2, 0.3, 0.4, 0.5 },
                                                  { 0.1, 0.2, 0.3, 0.4, 0.5 },
                                                  { 0.1, 0.2, 0.3, 0.4, 0.5 } };

  state.UnsafelySetCustomRateParameters(parameters);
  for (size_t i = 0; i < num_grid_cells; i++)
  {
    EXPECT_EQ(state.custom_rate_parameters_[i][0], 0.1);
    EXPECT_EQ(state.custom_rate_parameters_[i][1], 0.2);
    EXPECT_EQ(state.custom_rate_parameters_[i][2], 0.3);
    EXPECT_EQ(state.custom_rate_parameters_[i][3], 0.4);
    EXPECT_EQ(state.custom_rate_parameters_[i][4], 0.5);
  }
}

TEST(State, UnsafelySetCustomRateParameterCatchesTooFewGridCells)
{
  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = 2,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
  } };

  std::vector<std::vector<double>> parameters = { { 0.1, 0.2, 0.3, 0.4, 0.5 } };

  EXPECT_ANY_THROW(state.UnsafelySetCustomRateParameters(parameters));
}

TEST(State, UnsafelySetCustomRateParameterCatchesTooParameters)
{
  micm::State state{ micm::StateParameters{
      .number_of_grid_cells_ = 2,
      .number_of_rate_constants_ = 10,
      .variable_names_{ "foo", "bar", "baz", "quz" },
      .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
  } };

  std::vector<std::vector<double>> parameters = { { 0.1, 0.2, 0.3 } };

  EXPECT_ANY_THROW(state.UnsafelySetCustomRateParameters(parameters));
}