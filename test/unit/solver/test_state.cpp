#include <micm/solver/state.hpp>

#include <gtest/gtest.h>

TEST(State, DefaultConstructor)
{
  EXPECT_NO_THROW(micm::State state);
  EXPECT_NO_THROW(std::vector<micm::State<>> states(10));
}

TEST(State, Constructor)
{
  micm::State state{ micm::StateParameters{
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "quux", "corge" },
                     },
                     3 };

  EXPECT_EQ(state.conditions_.size(), 3);
  EXPECT_EQ(state.variable_map_["foo"], 0);
  EXPECT_EQ(state.variable_map_["bar"], 1);
  EXPECT_EQ(state.variable_map_["baz"], 2);
  EXPECT_EQ(state.variable_map_["quz"], 3);
  EXPECT_EQ(state.variables_.NumRows(), 3);
  EXPECT_EQ(state.variables_.NumColumns(), 4);
  EXPECT_EQ(state.variables_[0].Size(), 4);
  EXPECT_EQ(state.custom_rate_parameter_map_["quux"], 0);
  EXPECT_EQ(state.custom_rate_parameter_map_["corge"], 1);
  EXPECT_EQ(state.custom_rate_parameters_.NumRows(), 3);
  EXPECT_EQ(state.custom_rate_parameters_.NumColumns(), 2);
  EXPECT_EQ(state.custom_rate_parameters_[0].Size(), 2);
  EXPECT_EQ(state.rate_constants_.NumRows(), 3);
  EXPECT_EQ(state.rate_constants_.NumColumns(), 10);
  EXPECT_EQ(state.rate_constants_[0].Size(), 10);
}

TEST(State, CopyConstructor)
{
  micm::State original{ micm::StateParameters{ .number_of_rate_constants_ = 10,
                                               .variable_names_{ "foo", "bar", "baz", "quz" },
                                               .custom_rate_parameter_labels_{ "quux", "corge" },
                                               .relative_tolerance_ = 1e-05,
                                               .absolute_tolerance_ = { 1e-10, 1e-10, 1e-10, 1e-10 } },
                        3 };

  micm::State copy = original;

  EXPECT_EQ(copy.number_of_grid_cells_, original.number_of_grid_cells_);
  EXPECT_EQ(copy.variable_names_, original.variable_names_);
  EXPECT_EQ(copy.custom_rate_parameter_map_, original.custom_rate_parameter_map_);
  EXPECT_EQ(copy.variables_.NumRows(), original.variables_.NumRows());
  EXPECT_EQ(copy.variables_.NumColumns(), original.variables_.NumColumns());
  EXPECT_EQ(copy.custom_rate_parameters_.NumRows(), original.custom_rate_parameters_.NumRows());
  EXPECT_EQ(copy.custom_rate_parameters_.NumColumns(), original.custom_rate_parameters_.NumColumns());
  EXPECT_EQ(copy.rate_constants_.NumRows(), original.rate_constants_.NumRows());
  EXPECT_EQ(copy.rate_constants_.NumColumns(), original.rate_constants_.NumColumns());
  EXPECT_EQ(copy.relative_tolerance_, original.relative_tolerance_);
  EXPECT_EQ(copy.absolute_tolerance_, original.absolute_tolerance_);

  for (std::size_t i = 0; i < original.variables_.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < original.variables_.NumColumns(); ++j)
    {
      EXPECT_EQ(copy.variables_[i][j], original.variables_[i][j]);
    }
  }

  for (std::size_t i = 0; i < original.custom_rate_parameters_.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < original.custom_rate_parameters_.NumColumns(); ++j)
    {
      EXPECT_EQ(copy.custom_rate_parameters_[i][j], original.custom_rate_parameters_[i][j]);
    }
  }

  for (std::size_t i = 0; i < original.rate_constants_.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < original.rate_constants_.NumColumns(); ++j)
    {
      EXPECT_EQ(copy.rate_constants_[i][j], original.rate_constants_[i][j]);
    }
  }
}

TEST(State, CopyAssignmentOperator)
{
  micm::State original{ micm::StateParameters{ .number_of_rate_constants_ = 10,
                                               .variable_names_{ "foo", "bar", "baz", "quz" },
                                               .custom_rate_parameter_labels_{ "quux", "corge" },
                                               .relative_tolerance_ = 1e-05,
                                               .absolute_tolerance_ = { 1e-10, 1e-10, 1e-10, 1e-10 } },
                        3 };

  micm::State copy;
  copy = original;

  EXPECT_EQ(copy.number_of_grid_cells_, original.number_of_grid_cells_);
  EXPECT_EQ(copy.variable_names_, original.variable_names_);
  EXPECT_EQ(copy.custom_rate_parameter_map_, original.custom_rate_parameter_map_);
  EXPECT_EQ(copy.variables_.NumRows(), original.variables_.NumRows());
  EXPECT_EQ(copy.variables_.NumColumns(), original.variables_.NumColumns());
  EXPECT_EQ(copy.custom_rate_parameters_.NumRows(), original.custom_rate_parameters_.NumRows());
  EXPECT_EQ(copy.custom_rate_parameters_.NumColumns(), original.custom_rate_parameters_.NumColumns());
  EXPECT_EQ(copy.rate_constants_.NumRows(), original.rate_constants_.NumRows());
  EXPECT_EQ(copy.rate_constants_.NumColumns(), original.rate_constants_.NumColumns());
  EXPECT_EQ(copy.relative_tolerance_, original.relative_tolerance_);
  EXPECT_EQ(copy.absolute_tolerance_, original.absolute_tolerance_);

  for (std::size_t i = 0; i < original.variables_.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < original.variables_.NumColumns(); ++j)
    {
      EXPECT_EQ(copy.variables_[i][j], original.variables_[i][j]);
    }
  }

  for (std::size_t i = 0; i < original.custom_rate_parameters_.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < original.custom_rate_parameters_.NumColumns(); ++j)
    {
      EXPECT_EQ(copy.custom_rate_parameters_[i][j], original.custom_rate_parameters_[i][j]);
    }
  }

  for (std::size_t i = 0; i < original.rate_constants_.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < original.rate_constants_.NumColumns(); ++j)
    {
      EXPECT_EQ(copy.rate_constants_[i][j], original.rate_constants_[i][j]);
    }
  }
}

TEST(State, MoveConstructor)
{
  micm::State original{ micm::StateParameters{ .number_of_rate_constants_ = 10,
                                               .variable_names_{ "foo", "bar", "baz", "quz" },
                                               .custom_rate_parameter_labels_{ "quux", "corge" },
                                               .relative_tolerance_ = 1e-05,
                                               .absolute_tolerance_ = { 1e-10, 1e-10, 1e-10, 1e-10 } },
                        3 };

  auto expected_variable_names = original.variable_names_;
  auto expected_custom_rate_parameter_map = original.custom_rate_parameter_map_;
  auto expected_variables = original.variables_;
  auto expected_custom_rate_parameters = original.custom_rate_parameters_;
  auto expected_rate_constants = original.rate_constants_;

  micm::State moved = std::move(original);

  EXPECT_EQ(moved.number_of_grid_cells_, 3);
  EXPECT_EQ(moved.variable_names_, expected_variable_names);
  EXPECT_EQ(moved.custom_rate_parameter_map_, expected_custom_rate_parameter_map);
  EXPECT_EQ(moved.variables_.NumRows(), 3);
  EXPECT_EQ(moved.variables_.NumColumns(), 4);
  EXPECT_EQ(moved.custom_rate_parameters_.NumRows(), 3);
  EXPECT_EQ(moved.custom_rate_parameters_.NumColumns(), 2);
  EXPECT_EQ(moved.rate_constants_.NumRows(), 3);
  EXPECT_EQ(moved.rate_constants_.NumColumns(), 10);
  EXPECT_EQ(moved.relative_tolerance_, 1e-05);
  EXPECT_EQ(moved.absolute_tolerance_, std::vector<double>({ 1e-10, 1e-10, 1e-10, 1e-10 }));

  for (std::size_t i = 0; i < expected_variables.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < expected_variables.NumColumns(); ++j)
    {
      EXPECT_EQ(moved.variables_[i][j], expected_variables[i][j]);
    }
  }

  for (std::size_t i = 0; i < expected_custom_rate_parameters.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < expected_custom_rate_parameters.NumColumns(); ++j)
    {
      EXPECT_EQ(moved.custom_rate_parameters_[i][j], expected_custom_rate_parameters[i][j]);
    }
  }

  for (std::size_t i = 0; i < expected_rate_constants.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < expected_rate_constants.NumColumns(); ++j)
    {
      EXPECT_EQ(moved.rate_constants_[i][j], expected_rate_constants[i][j]);
    }
  }
}

TEST(State, MoveAssignmentOperator)
{
  micm::State original{ micm::StateParameters{ .number_of_rate_constants_ = 10,
                                               .variable_names_{ "foo", "bar", "baz", "quz" },
                                               .custom_rate_parameter_labels_{ "quux", "corge" },
                                               .relative_tolerance_ = 1e-05,
                                               .absolute_tolerance_ = { 1e-10, 1e-10, 1e-10, 1e-10 } },
                        3 };

  auto expected_variable_names = original.variable_names_;
  auto expected_custom_rate_parameter_map = original.custom_rate_parameter_map_;
  auto expected_variables = original.variables_;
  auto expected_custom_rate_parameters = original.custom_rate_parameters_;
  auto expected_rate_constants = original.rate_constants_;

  micm::State moved;
  moved = std::move(original);

  EXPECT_EQ(moved.number_of_grid_cells_, 3);
  EXPECT_EQ(moved.variable_names_, expected_variable_names);
  EXPECT_EQ(moved.custom_rate_parameter_map_, expected_custom_rate_parameter_map);
  EXPECT_EQ(moved.variables_.NumRows(), 3);
  EXPECT_EQ(moved.variables_.NumColumns(), 4);
  EXPECT_EQ(moved.custom_rate_parameters_.NumRows(), 3);
  EXPECT_EQ(moved.custom_rate_parameters_.NumColumns(), 2);
  EXPECT_EQ(moved.rate_constants_.NumRows(), 3);
  EXPECT_EQ(moved.rate_constants_.NumColumns(), 10);
  EXPECT_EQ(moved.relative_tolerance_, 1e-05);
  EXPECT_EQ(moved.absolute_tolerance_, std::vector<double>({ 1e-10, 1e-10, 1e-10, 1e-10 }));

  for (std::size_t i = 0; i < expected_variables.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < expected_variables.NumColumns(); ++j)
    {
      EXPECT_EQ(moved.variables_[i][j], expected_variables[i][j]);
    }
  }

  for (std::size_t i = 0; i < expected_custom_rate_parameters.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < expected_custom_rate_parameters.NumColumns(); ++j)
    {
      EXPECT_EQ(moved.custom_rate_parameters_[i][j], expected_custom_rate_parameters[i][j]);
    }
  }

  for (std::size_t i = 0; i < expected_rate_constants.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < expected_rate_constants.NumColumns(); ++j)
    {
      EXPECT_EQ(moved.rate_constants_[i][j], expected_rate_constants[i][j]);
    }
  }
}

TEST(State, SettingSingleConcentrationWithInvalidArgumentsThowsException)
{
  micm::State state{ micm::StateParameters{
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "quux", "corge" },
                     },
                     3 };
  EXPECT_ANY_THROW(state.SetConcentration(micm::Species{ "foo" }, 1.0));
  EXPECT_ANY_THROW(state.SetConcentration(micm::Species{ "foo" }, std::vector<double>{ 1.0, 2.0 }));
  EXPECT_ANY_THROW(state.SetConcentration(micm::Species{ "not foo" }, 1.0));
  EXPECT_ANY_THROW(state.SetConcentration(micm::Species{ "not foo" }, std::vector<double>{ 1.0, 2.0, 3.0, 4.0 }));
}

TEST(State, SetSingleConcentration)
{
  {
    micm::State state{ micm::StateParameters{
                           .number_of_rate_constants_ = 10,
                           .variable_names_{ "foo", "bar", "baz", "quz" },
                           .custom_rate_parameter_labels_{ "quux", "corge" },
                       },
                       3 };
    std::vector<double> concentrations{ 12.0, 42.0, 35.2 };
    state.SetConcentration(micm::Species{ "bar" }, concentrations);
    for (std::size_t i = 0; i < concentrations.size(); ++i)
      EXPECT_EQ(state.variables_[i][state.variable_map_["bar"]], concentrations[i]);
  }
  {
    micm::State state{ micm::StateParameters{
                           .number_of_rate_constants_ = 10,
                           .variable_names_{ "foo", "bar", "baz", "quz" },
                           .custom_rate_parameter_labels_{ "quux", "corge" },
                       },
                       1 };
    state.SetConcentration(micm::Species{ "bar" }, 324.2);
    EXPECT_EQ(state.variables_[0][state.variable_map_["bar"]], 324.2);
  }
}

TEST(State, SetConcentrationByElementSingleValue)
{
  std::string aitken_num_conc = "modal.aitken.number_concentration";
  std::string accum_num_conc = "modal.accumulation.number_concentration";

  micm::SystemParameters params;
  params.others_.push_back(aitken_num_conc);
  params.others_.push_back(accum_num_conc);

  micm::State state{ micm::StateParameters{
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz", aitken_num_conc, accum_num_conc },
                         .custom_rate_parameter_labels_{ "quux", "corge" },
                     },
                     1 };

  state.SetConcentration(aitken_num_conc, 42.0);
  state.SetConcentration(accum_num_conc, 12.0);

  EXPECT_EQ(state.variables_[0][state.variable_map_[aitken_num_conc]], 42.0);
  EXPECT_EQ(state.variables_[0][state.variable_map_[accum_num_conc]], 12.0);
}

TEST(State, SetConcentrationByElementVector)
{
  std::string aitken_num_conc = "modal.aitken.number_concentration";

  micm::SystemParameters params;
  params.others_.push_back(aitken_num_conc);

  micm::State state{ micm::StateParameters{
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz", aitken_num_conc },
                         .custom_rate_parameter_labels_{ "quux", "corge" },
                     },
                     3 };

  std::vector<double> concentrations{ 12.0, 42.0, 35.2 };

  state.SetConcentration(aitken_num_conc, concentrations);

  for (std::size_t i = 0; i < concentrations.size(); ++i)
    EXPECT_EQ(state.variables_[i][state.variable_map_[aitken_num_conc]], concentrations[i]);
}

TEST(State, SettingConcentrationsWithInvalidArguementsThrowsException)
{
  micm::State state{ micm::StateParameters{
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "quux", "corge" },
                     },
                     3 };

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
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "quux", "corge" },
                     },
                     num_grid_cells };

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
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
                     },
                     3 };

  // user input for custom rate parameters (unordered)
  std::unordered_map<std::string, std::vector<double>> custom_params = {
    { "O3", { 0.3 } }, { "O1", { 0.1 } }, { "O2", { 0.5 } }, { "BBB", { 0.7 } }, { "AAA", { 0.5 } }
  };

  EXPECT_ANY_THROW(state.SetCustomRateParameters(custom_params));
}

TEST(State, SettingCustomRateParameterWithInvalidLabelThrowsException)
{
  micm::State state{ micm::StateParameters{
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
                     },
                     1 };

  // user input for custom rate parameters (unordered)
  std::unordered_map<std::string, std::vector<double>> custom_params = {
    { "O3", { 0.3 } }, { "O1", { 0.1 } }, { "O2", { 0.5 } }, { "CCC", { 0.7 } }, { "AAA", { 0.5 } }
  };

  EXPECT_ANY_THROW(state.SetCustomRateParameters(custom_params));
}

TEST(State, SetCustomRateParameter)
{
  micm::State state{ micm::StateParameters{
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
                     },
                     1 };

  state.SetCustomRateParameter("O2", 42.3);
  EXPECT_EQ(state.custom_rate_parameters_[0][1], 42.3);
}

TEST(State, SetCustomRateParameters)
{
  uint32_t num_grid_cells = 3;

  micm::State state{ micm::StateParameters{
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
                     },
                     num_grid_cells };

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
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
                     },
                     1 };

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
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
                     },
                     num_grid_cells };

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
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
                     },
                     2 };

  std::vector<std::vector<double>> parameters = { { 0.1, 0.2, 0.3, 0.4, 0.5 } };

  EXPECT_ANY_THROW(state.UnsafelySetCustomRateParameters(parameters));
}

TEST(State, UnsafelySetCustomRateParameterCatchesTooParameters)
{
  micm::State state{ micm::StateParameters{
                         .number_of_rate_constants_ = 10,
                         .variable_names_{ "foo", "bar", "baz", "quz" },
                         .custom_rate_parameter_labels_{ "O1", "O2", "O3", "AAA", "BBB" },
                     },
                     2 };

  std::vector<std::vector<double>> parameters = { { 0.1, 0.2, 0.3 } };

  EXPECT_ANY_THROW(state.UnsafelySetCustomRateParameters(parameters));
}
