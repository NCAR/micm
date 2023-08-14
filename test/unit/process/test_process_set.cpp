#include <gtest/gtest.h>

#include <micm/process/process_set.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>
#include <random>

#include "test_process_set_policy.hpp"

template<class T>
using SparseMatrixTest = micm::SparseMatrix<T>;

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = micm::VectorMatrix<T, 4>;

template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;

TEST(ProcessSet, Matrix)
{
  testProcessSet<micm::Matrix, SparseMatrixTest, micm::ProcessSet>(
      [](const std::vector<micm::Process>& processes, const micm::State<micm::Matrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state };
      });
}

TEST(ProcessSet, VectorMatrix)
{
  testProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix, micm::ProcessSet>(
      [](const std::vector<micm::Process>& processes, const micm::State<Group1VectorMatrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state };
      });
  testProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix, micm::ProcessSet>(
      [](const std::vector<micm::Process>& processes, const micm::State<Group2VectorMatrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state };
      });
  testProcessSet<Group3VectorMatrix, Group3SparseVectorMatrix, micm::ProcessSet>(
      [](const std::vector<micm::Process>& processes, const micm::State<Group3VectorMatrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state };
      });
  testProcessSet<Group4VectorMatrix, Group4SparseVectorMatrix, micm::ProcessSet>(
      [](const std::vector<micm::Process>& processes, const micm::State<Group4VectorMatrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state };
      });
}

TEST(RandomProcessSet, Matrix)
{
  testRandomSystem<micm::Matrix, micm::ProcessSet>(
      2000,
      500,
      400,
      [](const std::vector<micm::Process>& processes, const micm::State<micm::Matrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state };
      });
  testRandomSystem<micm::Matrix, micm::ProcessSet>(
      3000,
      300,
      200,
      [](const std::vector<micm::Process>& processes, const micm::State<micm::Matrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state };
      });
  testRandomSystem<micm::Matrix, micm::ProcessSet>(
      4000,
      100,
      80,
      [](const std::vector<micm::Process>& processes, const micm::State<micm::Matrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state };
      });
}