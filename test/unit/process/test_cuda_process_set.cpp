#include <gtest/gtest.h>
#include <micm/process/process_set_cuda.cuh>
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

  size_t* number_of_reactants = set.number_of_reactants_.data();
  int number_of_reactants_size = set.number_of_reactants_.size();
  size_t* reactant_ids = set.reactant_ids_.data();
  int reactant_ids_size = set.reactant_ids_.size();
  size_t* number_of_products = set.number_of_products_.data();
  int number_of_products_size = set.number_of_products_.size();
  size_t* product_ids = set.product_ids_.data(); 
  int product_ids_size = set.product_ids_.size();
  double* yields = set.yields_.data();
  int yields_size = set.yields_.size(); 

  EXPECT_EQ(state.variables_.size(), 2);
  EXPECT_EQ(state.variables_[0].size(), 5);
  state.variables_[0] = { 0.1, 0.2, 0.3, 0.4, 0.5 };
  state.variables_[1] = { 1.1, 1.2, 1.3, 1.4, 1.5 };

  micm::Matrix<double> rate_constants{ 2, 3 };
  rate_constants[0] = { 10.0, 20.0, 30.0 };
  rate_constants[1] = { 110.0, 120.0, 130.0 };

  micm::Matrix<double> forcing{ 2, 5, 1000.0 };

  micm::cuda::AddForcingTerms_kernelSetup(
    number_of_reactants,
    number_of_reactants_size,
    reactant_ids,
    reactant_ids_size,
    number_of_products,
    number_of_products_size,
    product_ids,
    product_ids_size,
    yields,
    yields_size,
    rate_constants, 
    state.variables_, 
    forcing);

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

//   auto non_zero_elements = set.NonZeroJacobianElements();

//   // ---- foo  bar  baz  quz  quuz
//   // foo   0    1    2    -    -
//   // bar   3    4    5    -    -
//   // baz   6    -    7    -    -
//   // quz   -    8    -    9    -
//   // quuz 10    -   11    -    -

//   auto elem = non_zero_elements.begin();
//   compare_pair(*elem, index_pair(0, 0));
//   compare_pair(*(++elem), index_pair(0, 1));
//   compare_pair(*(++elem), index_pair(0, 2));
//   compare_pair(*(++elem), index_pair(1, 0));
//   compare_pair(*(++elem), index_pair(1, 1));
//   compare_pair(*(++elem), index_pair(1, 2));
//   compare_pair(*(++elem), index_pair(2, 0));
//   compare_pair(*(++elem), index_pair(2, 2));
//   compare_pair(*(++elem), index_pair(3, 1));
//   compare_pair(*(++elem), index_pair(3, 3));
//   compare_pair(*(++elem), index_pair(4, 0));
//   compare_pair(*(++elem), index_pair(4, 2));

//   auto builder = micm::SparseMatrix<double>::create(5).number_of_blocks(2).initial_value(100.0);
//   for (auto& elem : non_zero_elements)
//     builder = builder.with_element(elem.first, elem.second);
//   micm::SparseMatrix<double> jacobian{ builder };
//   set.SetJacobianFlatIds(jacobian);
//   set.AddJacobianTerms(rate_constants, state.variables_, jacobian);

//   EXPECT_EQ(jacobian[0][0][0], 100.0 - 10.0 * 0.3);  // foo -> foo
//   EXPECT_EQ(jacobian[1][0][0], 100.0 - 110.0 * 1.3);
//   EXPECT_EQ(jacobian[0][0][1], 100.0 + 20.0);  // foo -> bar
//   EXPECT_EQ(jacobian[1][0][1], 100.0 + 120.0);
//   EXPECT_EQ(jacobian[0][0][2], 100.0 - 10.0 * 0.1);  // foo -> baz
//   EXPECT_EQ(jacobian[1][0][2], 100.0 - 110.0 * 1.1);
//   EXPECT_EQ(jacobian[0][1][0], 100.0 + 10.0 * 0.3);  // bar -> foo
//   EXPECT_EQ(jacobian[1][1][0], 100.0 + 110.0 * 1.3);
//   EXPECT_EQ(jacobian[0][1][1], 100.0 - 20.0);  // bar -> bar
//   EXPECT_EQ(jacobian[1][1][1], 100.0 - 120.0);
//   EXPECT_EQ(jacobian[0][1][2], 100.0 + 10.0 * 0.1);  // bar -> baz
//   EXPECT_EQ(jacobian[1][1][2], 100.0 + 110.0 * 1.1);
//   EXPECT_EQ(jacobian[0][2][0], 100.0 - 10.0 * 0.3);  // baz -> foo
//   EXPECT_EQ(jacobian[1][2][0], 100.0 - 110.0 * 1.3);
//   EXPECT_EQ(jacobian[0][2][2], 100.0 - 10.0 * 0.1);  // baz -> baz
//   EXPECT_EQ(jacobian[1][2][2], 100.0 - 110.0 * 1.1);
//   EXPECT_EQ(jacobian[0][3][1], 100.0 + 1.4 * 20.0);  // quz -> bar
//   EXPECT_EQ(jacobian[1][3][1], 100.0 + 1.4 * 120.0);
//   EXPECT_EQ(jacobian[0][3][3], 100.0 - 30.0);  // quz -> quz
//   EXPECT_EQ(jacobian[1][3][3], 100.0 - 130.0);
//   EXPECT_EQ(jacobian[0][4][0], 100.0 + 2.4 * 10.0 * 0.3);  // quuz -> foo
//   EXPECT_EQ(jacobian[1][4][0], 100.0 + 2.4 * 110.0 * 1.3);
//   EXPECT_EQ(jacobian[0][4][2], 100.0 + 2.4 * 10.0 * 0.1);  // quuz -> baz
//   EXPECT_EQ(jacobian[1][4][2], 100.0 + 2.4 * 110.0 * 1.1);
 }








// // TEST(ProcessSet, Constructor)
// // {
// //   // micm::cuda::AddForcingTerms_kernelSetup(...);
// // }