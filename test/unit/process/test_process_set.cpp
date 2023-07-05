#include <gtest/gtest.h>

#include <micm/process/process.hpp>
#include <micm/process/process_set.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/vector_matrix.hpp>
#include <random>

using yields = std::pair<micm::Species, double>;
using index_pair = std::pair<std::size_t, std::size_t>;

void compare_pair(const index_pair& a, const index_pair& b)
{
  EXPECT_EQ(a.first, b.first);
  EXPECT_EQ(a.second, b.second);
}

template<template<class> class MatrixPolicy>
void testProcessSet()
{
  auto foo = micm::Species("foo");
  auto bar = micm::Species("bar");
  auto baz = micm::Species("baz");
  auto quz = micm::Species("quz");
  auto quuz = micm::Species("quuz");

  micm::Phase gas_phase{ std::vector<micm::Species>{ foo, bar, baz, quz, quuz } };

  micm::State<MatrixPolicy> state{ micm::StateParameters{ .state_variable_names_{ "foo", "bar", "baz", "quz", "quuz" },
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

  MatrixPolicy<double> rate_constants{ 2, 3 };
  rate_constants[0] = { 10.0, 20.0, 30.0 };
  rate_constants[1] = { 110.0, 120.0, 130.0 };

  MatrixPolicy<double> forcing{ 2, 5, 1000.0 };

  set.AddForcingTerms<MatrixPolicy>(rate_constants, state.variables_, forcing);

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
  // foo   0    1    2    -    -
  // bar   3    4    5    -    -
  // baz   6    -    7    -    -
  // quz   -    8    -    9    -
  // quuz 10    -   11    -    -

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

  auto builder = micm::SparseMatrix<double>::create(5).number_of_blocks(2).initial_value(100.0);
  for (auto& elem : non_zero_elements)
    builder = builder.with_element(elem.first, elem.second);
  micm::SparseMatrix<double> jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);
  set.AddJacobianTerms<MatrixPolicy>(rate_constants, state.variables_, jacobian);

  EXPECT_EQ(jacobian[0][0][0], 100.0 - 10.0 * 0.3);  // foo -> foo
  EXPECT_EQ(jacobian[1][0][0], 100.0 - 110.0 * 1.3);
  EXPECT_EQ(jacobian[0][0][1], 100.0 + 20.0);  // foo -> bar
  EXPECT_EQ(jacobian[1][0][1], 100.0 + 120.0);
  EXPECT_EQ(jacobian[0][0][2], 100.0 - 10.0 * 0.1);  // foo -> baz
  EXPECT_EQ(jacobian[1][0][2], 100.0 - 110.0 * 1.1);
  EXPECT_EQ(jacobian[0][1][0], 100.0 + 10.0 * 0.3);  // bar -> foo
  EXPECT_EQ(jacobian[1][1][0], 100.0 + 110.0 * 1.3);
  EXPECT_EQ(jacobian[0][1][1], 100.0 - 20.0);  // bar -> bar
  EXPECT_EQ(jacobian[1][1][1], 100.0 - 120.0);
  EXPECT_EQ(jacobian[0][1][2], 100.0 + 10.0 * 0.1);  // bar -> baz
  EXPECT_EQ(jacobian[1][1][2], 100.0 + 110.0 * 1.1);
  EXPECT_EQ(jacobian[0][2][0], 100.0 - 10.0 * 0.3);  // baz -> foo
  EXPECT_EQ(jacobian[1][2][0], 100.0 - 110.0 * 1.3);
  EXPECT_EQ(jacobian[0][2][2], 100.0 - 10.0 * 0.1);  // baz -> baz
  EXPECT_EQ(jacobian[1][2][2], 100.0 - 110.0 * 1.1);
  EXPECT_EQ(jacobian[0][3][1], 100.0 + 1.4 * 20.0);  // quz -> bar
  EXPECT_EQ(jacobian[1][3][1], 100.0 + 1.4 * 120.0);
  EXPECT_EQ(jacobian[0][3][3], 100.0 - 30.0);  // quz -> quz
  EXPECT_EQ(jacobian[1][3][3], 100.0 - 130.0);
  EXPECT_EQ(jacobian[0][4][0], 100.0 + 2.4 * 10.0 * 0.3);  // quuz -> foo
  EXPECT_EQ(jacobian[1][4][0], 100.0 + 2.4 * 110.0 * 1.3);
  EXPECT_EQ(jacobian[0][4][2], 100.0 + 2.4 * 10.0 * 0.1);  // quuz -> baz
  EXPECT_EQ(jacobian[1][4][2], 100.0 + 2.4 * 110.0 * 1.1);
}

TEST(ProcessSet, Matrix)
{
  testProcessSet<micm::Matrix>();
}

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = micm::VectorMatrix<T, 4>;

TEST(ProcessSet, VectorMatrix)
{
  testProcessSet<Group1VectorMatrix>();
  testProcessSet<Group2VectorMatrix>();
  testProcessSet<Group3VectorMatrix>();
  testProcessSet<Group4VectorMatrix>();
}

template<template<class> class MatrixPolicy>
void testRandomSystem(std::size_t n_cells, std::size_t n_reactions, std::size_t n_species)
{
  auto get_n_react = std::bind(std::uniform_int_distribution<>(0, 3), std::default_random_engine());
  auto get_n_product = std::bind(std::uniform_int_distribution<>(0, 10), std::default_random_engine());
  auto get_species_id = std::bind(std::uniform_int_distribution<>(0, n_species - 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 4.0), std::default_random_engine());

  std::vector<micm::Species> species{};
  std::vector<std::string> species_names{};
  for (std::size_t i = 0; i < n_species; ++i)
  {
    species.push_back(micm::Species{ std::to_string(i) });
    species_names.push_back(std::to_string(i));
  }
  micm::Phase gas_phase{ species };
  micm::State<MatrixPolicy> state{ micm::StateParameters{ .state_variable_names_{ species_names },
                                                          .number_of_grid_cells_ = n_cells,
                                                          .number_of_custom_parameters_ = 0,
                                                          .number_of_rate_constants_ = n_reactions } };

  std::vector<micm::Process> processes{};
  for (std::size_t i = 0; i < n_reactions; ++i)
  {
    auto n_react = get_n_react();
    std::vector<micm::Species> reactants{};
    for (std::size_t i_react = 0; i_react < n_react; ++i_react)
    {
      reactants.push_back({ std::to_string(get_species_id()) });
    }
    auto n_product = get_n_product();
    std::vector<yields> products{};
    for (std::size_t i_prod = 0; i_prod < n_product; ++i_prod)
    {
      products.push_back(yields(std::to_string(get_species_id()), 1.2));
    }
    processes.push_back(micm::Process::create().reactants(reactants).products(products).phase(gas_phase));
  }

  micm::ProcessSet set{ processes, state };

  for (auto& elem : state.variables_.AsVector())
    elem = get_double();

  MatrixPolicy<double> rate_constants{ n_cells, n_reactions };
  for (auto& elem : rate_constants.AsVector())
    elem = get_double();

  MatrixPolicy<double> forcing{ n_cells, n_species, 1000.0 };

  set.AddForcingTerms<MatrixPolicy>(rate_constants, state.variables_, forcing);
}

TEST(RandomProcessSet, Matrix)
{
  testRandomSystem<micm::Matrix>(2000, 500, 400);
  testRandomSystem<micm::Matrix>(3000, 300, 200);
  testRandomSystem<micm::Matrix>(4000, 100, 80);
}