#include <gtest/gtest.h>

#include <micm/process/process.hpp>
#include <micm/process/process_set.hpp>

using yields = std::pair<micm::Species, double>;
using index_pair = std::pair<std::size_t, std::size_t>;

void compare_pair(const index_pair& a, const index_pair& b)
{
  EXPECT_EQ(a.first, b.first);
  EXPECT_EQ(a.second, b.second);
}

TEST(ProcessSet, Constructor)
{
  auto foo = micm::Species("foo");
  auto bar = micm::Species("bar");
  auto baz = micm::Species("baz");
  auto quz = micm::Species("quz");
  auto quuz = micm::Species("quuz");

  micm::Phase gas_phase{ std::vector<micm::Species>{ foo, bar, baz, quz, quuz } };

  micm::State state{ micm::StateParameters{ .state_variable_names_{ "foo", "bar", "baz", "quz", "quuz" },
                                            .number_of_grid_cells_ = 2,
                                            .number_of_custom_parameters_ = 0,
                                            .number_of_rate_constants_ = 3 } };

  micm::Process r1 =
      micm::Process::create().reactants({ foo, baz }).products({ yields(bar, 1), yields(quuz, 2.4) }).phase(gas_phase);

  micm::Process r2 =
      micm::Process::create().reactants({ bar }).products({ yields(foo, 1), yields(quz, 1.4) }).phase(gas_phase);

  micm::Process r3 = micm::Process::create().reactants({ quz }).products({}).phase(gas_phase);

  micm::ProcessSet set{ std::vector<micm::Process>{ r1, r2, r3 }, state };

  EXPECT_EQ(state.variables_.size(), 2);
  EXPECT_EQ(state.variables_[0].size(), 5);
  state.variables_[0] = { 0.1, 0.2, 0.3, 0.4, 0.5 };
  state.variables_[1] = { 1.1, 1.2, 1.3, 1.4, 1.5 };

  micm::Matrix<double> rate_constants{ 2, 3 };
  rate_constants[0] = { 10.0, 20.0, 30.0 };
  rate_constants[1] = { 110.0, 120.0, 130.0 };

  micm::Matrix<double> forcing{ 2, 5, 1000.0 };

  set.AddForcingTerms(rate_constants, state.variables_, forcing);

  EXPECT_EQ(forcing[0][0], 1000.0 - 10.0 * 0.1 * 0.3 + 20.0 * 0.2);
  EXPECT_EQ(forcing[1][0], 1000.0 - 110.0 * 1.1 * 1.3 + 120.0 * 1.2);
  EXPECT_EQ(forcing[0][1], 1000.0 + 10.0 * 0.1 * 0.3 - 20.0 * 0.2);
  EXPECT_EQ(forcing[1][1], 1000.0 + 110.0 * 1.1 * 1.3 - 120.0 * 1.2);
  EXPECT_EQ(forcing[0][2], 1000.0 - 10.0 * 0.1 * 0.3);
  EXPECT_EQ(forcing[1][2], 1000.0 - 110.0 * 1.1 * 1.3);
  EXPECT_EQ(forcing[0][3], 1000.0 + 20.0 * 0.2 * 1.4 - 30.0 * 0.4);
  EXPECT_EQ(forcing[1][3], 1000.0 + 120.0 * 1.2 * 1.4 - 130.0 * 1.4);
  EXPECT_EQ(forcing[0][4], 1000.0 + 10.0 * 0.1 * 0.3 * 2.4);
  EXPECT_EQ(forcing[1][4], 1000.0 + 110.0 * 1.1 * 1.3 * 2.4);

  auto non_zero_elements = set.NonZeroJacobianElements();

  // ---- foo  bar  baz  quz  quuz
  // foo   X    X    X    0    0
  // bar   X    X    X    0    0
  // baz   X    0    X    0    0
  // quz   0    X    0    X    0
  // quuz  X    0    X    0    0

  auto elem = non_zero_elements.begin();
  compare_pair(*elem, index_pair(0, 0));
  compare_pair(*(++elem), index_pair(0, 1));
  compare_pair(*(++elem), index_pair(0, 2));
  compare_pair(*(++elem), index_pair(1, 0));
  compare_pair(*(++elem), index_pair(1, 1));
  compare_pair(*(++elem), index_pair(1, 2));
  compare_pair(*(++elem), index_pair(2, 0));
  compare_pair(*(++elem), index_pair(2, 2));
  compare_pair(*(++elem), index_pair(3, 1));
  compare_pair(*(++elem), index_pair(3, 3));
  compare_pair(*(++elem), index_pair(4, 0));
  compare_pair(*(++elem), index_pair(4, 2));
}