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
      [](const std::vector<micm::Process>& processes,
         const micm::State<micm::Matrix, SparseMatrixTest>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state.variable_map_ };
      });
}

TEST(ProcessSet, VectorMatrix)
{
  testProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix, micm::ProcessSet>(
      [](const std::vector<micm::Process>& processes,
         const micm::State<Group1VectorMatrix, Group1SparseVectorMatrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state.variable_map_ };
      });
  testProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix, micm::ProcessSet>(
      [](const std::vector<micm::Process>& processes,
         const micm::State<Group2VectorMatrix, Group2SparseVectorMatrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state.variable_map_ };
      });
  testProcessSet<Group3VectorMatrix, Group3SparseVectorMatrix, micm::ProcessSet>(
      [](const std::vector<micm::Process>& processes,
         const micm::State<Group3VectorMatrix, Group3SparseVectorMatrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state.variable_map_ };
      });
  testProcessSet<Group4VectorMatrix, Group4SparseVectorMatrix, micm::ProcessSet>(
      [](const std::vector<micm::Process>& processes,
         const micm::State<Group4VectorMatrix, Group4SparseVectorMatrix>& state) -> micm::ProcessSet {
        return micm::ProcessSet{ processes, state.variable_map_ };
      });
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