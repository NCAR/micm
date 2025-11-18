#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/process.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

#include <random>

using namespace micm;
using index_pair = std::pair<std::size_t, std::size_t>;

void compare_pair(const index_pair& a, const index_pair& b)
{
  EXPECT_EQ(a.first, b.first);
  EXPECT_EQ(a.second, b.second);
}

template<class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy>
void testProcessSet()
{
  auto foo = Species("foo");
  auto bar = Species("bar");
  auto baz = Species("baz");
  auto quz = Species("quz");
  auto quuz = Species("quuz");
  auto qux = Species("qux");
  auto corge = Species("corge");
  qux.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.72; };

  Phase gas_phase{ "gas", std::vector<PhaseSpecies>{ foo, bar, baz, quz, quuz, corge } };
  State<DenseMatrixPolicy, SparseMatrixPolicy> state(
      StateParameters{ .number_of_rate_constants_ = 3, .variable_names_{ "foo", "bar", "baz", "quz", "quuz", "corge" } }, 2);

  ArrheniusRateConstant arrhenius_rate_constant({ .A_ = 12.2, .C_ = 300.0 });

  Process r1 = ChemicalReactionBuilder()
                   .SetReactants({ foo, baz })
                   .SetProducts({ Yield(bar, 1), Yield(quuz, 2.4) })
                   .SetRateConstant(arrhenius_rate_constant)
                   .SetPhase(gas_phase)
                   .Build();

  Process r2 = ChemicalReactionBuilder()
                   .SetReactants({ bar, qux })
                   .SetProducts({ Yield(foo, 1), Yield(quz, 1.4) })
                   .SetRateConstant(arrhenius_rate_constant)
                   .SetPhase(gas_phase)
                   .Build();

  Process r3 = ChemicalReactionBuilder()
                   .SetReactants({ quz })
                   .SetProducts({})
                   .SetRateConstant(arrhenius_rate_constant)
                   .SetPhase(gas_phase)
                   .Build();

  Process r4 = ChemicalReactionBuilder()
                   .SetReactants({ baz, qux })
                   .SetProducts({ Yield(bar, 1), Yield(quz, 2.5) })
                   .SetRateConstant(arrhenius_rate_constant)
                   .SetPhase(gas_phase)
                   .Build();

  auto used_species = RatesPolicy::SpeciesUsed(std::vector<Process>{ r1, r2, r3, r4 });

  EXPECT_EQ(used_species.size(), 6);
  EXPECT_TRUE(used_species.contains("foo"));
  EXPECT_TRUE(used_species.contains("bar"));
  EXPECT_TRUE(used_species.contains("baz"));
  EXPECT_TRUE(used_species.contains("quz"));
  EXPECT_TRUE(used_species.contains("quuz"));
  EXPECT_TRUE(used_species.contains("qux"));
  EXPECT_FALSE(used_species.contains("corge"));

  RatesPolicy set = RatesPolicy(std::vector<Process>{ r1, r2, r3, r4 }, state.variable_map_);

  EXPECT_EQ(state.variables_.NumRows(), 2);
  EXPECT_EQ(state.variables_.NumColumns(), 6);
  state.variables_[0] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.0 };
  state.variables_[1] = { 1.1, 1.2, 1.3, 1.4, 1.5, 0.0 };
  state.conditions_[0].air_density_ = 70.0;
  state.conditions_[1].air_density_ = 80.0;
  DenseMatrixPolicy rate_constants{ 2, 4 };
  // the rate constants will have been calculated and combined with
  // parameterized species before calculating forcing terms
  rate_constants[0] = { 10.0, 20.0 * 70.0 * 0.72, 30.0, 40.0 * 70.0 * 0.72 };
  rate_constants[1] = { 110.0, 120.0 * 80.0 * 0.72, 130.0, 140.0 * 80.0 * 0.72 };

  // Copy input-only variables to the device
  CheckCopyToDevice<DenseMatrixPolicy>(rate_constants);
  CheckCopyToDevice<DenseMatrixPolicy>(state.variables_);

  DenseMatrixPolicy forcing{ 2, 5, 1000.0 };

  CheckCopyToDevice<DenseMatrixPolicy>(forcing);

  set.template AddForcingTerms<DenseMatrixPolicy>(rate_constants, state.variables_, forcing);

  CheckCopyToHost<DenseMatrixPolicy>(forcing);

  EXPECT_DOUBLE_EQ(forcing[0][0], 1000.0 - 10.0 * 0.1 * 0.3 + 20.0 * 70.0 * 0.72 * 0.2);  // foo
  EXPECT_DOUBLE_EQ(forcing[1][0], 1000.0 - 110.0 * 1.1 * 1.3 + 120.0 * 80.0 * 0.72 * 1.2);
  EXPECT_DOUBLE_EQ(forcing[0][1], 1000.0 + 10.0 * 0.1 * 0.3 - 20.0 * 0.2 * 70.0 * 0.72 + 40.0 * 70.0 * 0.72 * 0.3);  // bar
  EXPECT_DOUBLE_EQ(forcing[1][1], 1000.0 + 110.0 * 1.1 * 1.3 - 120.0 * 1.2 * 80.0 * 0.72 + 140.0 * 80.0 * 0.72 * 1.3);
  EXPECT_DOUBLE_EQ(forcing[0][2], 1000.0 - 10.0 * 0.1 * 0.3 - 40.0 * 70.0 * 0.72 * 0.3);  // baz
  EXPECT_DOUBLE_EQ(forcing[1][2], 1000.0 - 110.0 * 1.1 * 1.3 - 140.0 * 80.0 * 0.72 * 1.3);
  EXPECT_DOUBLE_EQ(
      forcing[0][3], 1000.0 + 20.0 * 70.0 * 0.72 * 0.2 * 1.4 - 30.0 * 0.4 + 40.0 * 70.0 * 0.72 * 2.5 * 0.3);  // quz
  EXPECT_DOUBLE_EQ(forcing[1][3], 1000.0 + 120.0 * 80.0 * 0.72 * 1.2 * 1.4 - 130.0 * 1.4 + 140.0 * 80.0 * 0.72 * 2.5 * 1.3);
  EXPECT_DOUBLE_EQ(forcing[0][4], 1000.0 + 10.0 * 0.1 * 0.3 * 2.4);  // quuz
  EXPECT_DOUBLE_EQ(forcing[1][4], 1000.0 + 110.0 * 1.1 * 1.3 * 2.4);

  auto non_zero_elements = set.NonZeroJacobianElements();
  // ---- foo  bar  baz  quz  quuz
  // foo   0    1    2    -    -
  // bar   3    4    5    -    -
  // baz   6    -    7    -    -
  // quz   -    8    9   10    -
  // quuz 11    -   12    -    -

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
  compare_pair(*(++elem), index_pair(3, 2));
  compare_pair(*(++elem), index_pair(3, 3));
  compare_pair(*(++elem), index_pair(4, 0));
  compare_pair(*(++elem), index_pair(4, 2));

  auto builder = SparseMatrixPolicy::Create(5).SetNumberOfBlocks(2).InitialValue(100.0);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);
  SparseMatrixPolicy jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  CheckCopyToDevice<SparseMatrixPolicy>(jacobian);

  set.SubtractJacobianTerms(rate_constants, state.variables_, jacobian);

  CheckCopyToHost<SparseMatrixPolicy>(jacobian);

  EXPECT_DOUBLE_EQ(jacobian[0][0][0], 100.0 + 10.0 * 0.3);  // foo -> foo
  EXPECT_DOUBLE_EQ(jacobian[1][0][0], 100.0 + 110.0 * 1.3);
  EXPECT_DOUBLE_EQ(jacobian[0][0][1], 100.0 - 20.0 * 70.0 * 0.72);  // foo -> bar
  EXPECT_DOUBLE_EQ(jacobian[1][0][1], 100.0 - 120.0 * 80.0 * 0.72);
  EXPECT_DOUBLE_EQ(jacobian[0][0][2], 100.0 + 10.0 * 0.1);  // foo -> baz
  EXPECT_DOUBLE_EQ(jacobian[1][0][2], 100.0 + 110.0 * 1.1);
  EXPECT_DOUBLE_EQ(jacobian[0][1][0], 100.0 - 10.0 * 0.3);  // bar -> foo
  EXPECT_DOUBLE_EQ(jacobian[1][1][0], 100.0 - 110.0 * 1.3);
  EXPECT_DOUBLE_EQ(jacobian[0][1][1], 100.0 + 20.0 * 70.0 * 0.72);  // bar -> bar
  EXPECT_DOUBLE_EQ(jacobian[1][1][1], 100.0 + 120.0 * 80.0 * 0.72);
  EXPECT_DOUBLE_EQ(jacobian[0][1][2], 100.0 - 10.0 * 0.1 - 40.0 * 70.0 * 0.72);  // bar -> baz
  EXPECT_DOUBLE_EQ(jacobian[1][1][2], 100.0 - 110.0 * 1.1 - 140.0 * 80.0 * 0.72);
  EXPECT_DOUBLE_EQ(jacobian[0][2][0], 100.0 + 10.0 * 0.3);  // baz -> foo
  EXPECT_DOUBLE_EQ(jacobian[1][2][0], 100.0 + 110.0 * 1.3);
  EXPECT_DOUBLE_EQ(jacobian[0][2][2], 100.0 + 10.0 * 0.1 + 40.0 * 70.0 * 0.72);  // baz -> baz
  EXPECT_DOUBLE_EQ(jacobian[1][2][2], 100.0 + 110.0 * 1.1 + 140.0 * 80.0 * 0.72);
  EXPECT_DOUBLE_EQ(jacobian[0][3][1], 100.0 - 1.4 * 20.0 * 70.0 * 0.72);  // quz -> bar
  EXPECT_DOUBLE_EQ(jacobian[1][3][1], 100.0 - 1.4 * 120.0 * 80.0 * 0.72);
  EXPECT_DOUBLE_EQ(jacobian[0][3][2], 100.0 - 2.5 * 40.0 * 70.0 * 0.72);  // quz -> baz
  EXPECT_DOUBLE_EQ(jacobian[1][3][2], 100.0 - 2.5 * 140.0 * 80.0 * 0.72);
  EXPECT_DOUBLE_EQ(jacobian[0][3][3], 100.0 + 30.0);  // quz -> quz
  EXPECT_DOUBLE_EQ(jacobian[1][3][3], 100.0 + 130.0);
  EXPECT_DOUBLE_EQ(jacobian[0][4][0], 100.0 - 2.4 * 10.0 * 0.3);  // quuz -> foo
  EXPECT_DOUBLE_EQ(jacobian[1][4][0], 100.0 - 2.4 * 110.0 * 1.3);
  EXPECT_DOUBLE_EQ(jacobian[0][4][2], 100.0 - 2.4 * 10.0 * 0.1);  // quuz -> baz
  EXPECT_DOUBLE_EQ(jacobian[1][4][2], 100.0 - 2.4 * 110.0 * 1.1);
}

template<class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy>
void testRandomSystem(std::size_t n_cells, std::size_t n_reactions, std::size_t n_species)
{
  auto get_n_react = std::bind(std::uniform_int_distribution<>(0, 3), std::default_random_engine());
  auto get_n_product = std::bind(std::uniform_int_distribution<>(0, 10), std::default_random_engine());
  auto get_species_id = std::bind(std::uniform_int_distribution<>(0, n_species - 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 4.0), std::default_random_engine());

  std::vector<PhaseSpecies> phase_species{};
  std::vector<std::string> species_names{};
  phase_species.reserve(n_species);
  species_names.reserve(n_species);
  for (std::size_t i = 0; i < n_species; ++i)
  {
    phase_species.emplace_back(PhaseSpecies(Species(std::to_string(i))));
    species_names.emplace_back(std::to_string(i));
  }
  Phase gas_phase{ "gas", phase_species };

  ArrheniusRateConstant arrhenius_rate_constant({ .A_ = 12.2, .C_ = 300.0 });
  State<DenseMatrixPolicy, SparseMatrixPolicy> state{ StateParameters{
                                                          .number_of_rate_constants_ = n_reactions,
                                                          .variable_names_{ species_names },
                                                      },
                                                      n_cells };
  std::vector<Process> processes{};
  for (std::size_t i = 0; i < n_reactions; ++i)
  {
    auto n_react = get_n_react();
    std::vector<Species> reactants{};
    for (std::size_t i_react = 0; i_react < n_react; ++i_react)
    {
      reactants.push_back({ std::to_string(get_species_id()) });
    }
    auto n_product = get_n_product();
    std::vector<Yield> products{};
    for (std::size_t i_prod = 0; i_prod < n_product; ++i_prod)
    {
      products.push_back(Yield(std::to_string(get_species_id()), 1.2));
    }
    auto proc = ChemicalReactionBuilder()
                    .SetReactants(reactants)
                    .SetProducts(products)
                    .SetRateConstant(arrhenius_rate_constant)
                    .SetPhase(gas_phase)
                    .Build();
    processes.push_back(proc);
  }
  RatesPolicy set = RatesPolicy(processes, state.variable_map_);

  for (auto& elem : state.variables_.AsVector())
    elem = get_double();

  DenseMatrixPolicy rate_constants{ n_cells, n_reactions };
  for (auto& elem : rate_constants.AsVector())
    elem = get_double();
  DenseMatrixPolicy forcing{ n_cells, n_species, 1000.0 };

  CheckCopyToDevice<DenseMatrixPolicy>(forcing);

  set.template AddForcingTerms<DenseMatrixPolicy>(rate_constants, state.variables_, forcing);

  CheckCopyToHost<DenseMatrixPolicy>(forcing);
}
