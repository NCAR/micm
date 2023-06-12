#include <gtest/gtest.h>

#include <micm/solver/state.hpp>

TEST(State, Constructor)
{
  micm::StateParameters params;
  params.state_variable_names_ = { "foo", "bar", "baz", "quz" };
  params.number_of_grid_cells_ = 3;
  params.number_of_custom_parameters_ = 5;
  params.number_of_rate_constants_ = 10;

  micm::State state(params);


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