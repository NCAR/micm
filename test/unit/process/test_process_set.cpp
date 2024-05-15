#include "test_process_set_policy.hpp"

#include <micm/process/process_set.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <random>

using SparseMatrixTest = micm::SparseMatrix<double>;

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = micm::VectorMatrix<T, 4>;

using Group1SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<4>>;

TEST(ProcessSet, Matrix)
{
  testProcessSet<micm::Matrix, SparseMatrixTest, micm::ProcessSet>();
}

TEST(ProcessSet, VectorMatrix)
{
  testProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix, micm::ProcessSet>();
  testProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix, micm::ProcessSet>();
  testProcessSet<Group3VectorMatrix, Group3SparseVectorMatrix, micm::ProcessSet>();
  testProcessSet<Group4VectorMatrix, Group4SparseVectorMatrix, micm::ProcessSet>();
}

TEST(RandomProcessSet, Matrix)
{
  testRandomSystem<micm::Matrix, SparseMatrixTest, micm::ProcessSet>(
      200,
      50,
      40,
      [](const std::vector<micm::Process>& processes,
         const micm::State<micm::Matrix, SparseMatrixTest>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state.variable_map_ };
      });
  testRandomSystem<micm::Matrix, SparseMatrixTest, micm::ProcessSet>(
      300,
      30,
      20,
      [](const std::vector<micm::Process>& processes,
         const micm::State<micm::Matrix, SparseMatrixTest>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state.variable_map_ };
      });
  testRandomSystem<micm::Matrix, SparseMatrixTest, micm::ProcessSet>(
      400,
      100,
      80,
      [](const std::vector<micm::Process>& processes,
         const micm::State<micm::Matrix, SparseMatrixTest>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state.variable_map_ };
      });
}