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

using Group1VectorMatrix = micm::VectorMatrix<double, 1>;
using Group2VectorMatrix = micm::VectorMatrix<double, 2>;
using Group3VectorMatrix = micm::VectorMatrix<double, 3>;
using Group4VectorMatrix = micm::VectorMatrix<double, 4>;

using Group1SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<4>>;

TEST(ProcessSet, Matrix)
{
  testProcessSet<micm::Matrix<double>, SparseMatrixTest, micm::ProcessSet<micm::Matrix<double>, SparseMatrixTest>>();
}

TEST(ProcessSet, VectorMatrix)
{
  testProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix, micm::ProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  testProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix, micm::ProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  testProcessSet<Group3VectorMatrix, Group3SparseVectorMatrix, micm::ProcessSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  testProcessSet<Group4VectorMatrix, Group4SparseVectorMatrix, micm::ProcessSet<Group4VectorMatrix, Group4SparseVectorMatrix>>();
}

TEST(RandomProcessSet, Matrix)
{
  testRandomSystem<micm::Matrix<double>, SparseMatrixTest, micm::ProcessSet<micm::Matrix<double>, SparseMatrixTest>>(200, 50, 40);
  testRandomSystem<micm::Matrix<double>, SparseMatrixTest, micm::ProcessSet<micm::Matrix<double>, SparseMatrixTest>>(300, 30, 20);
  testRandomSystem<micm::Matrix<double>, SparseMatrixTest, micm::ProcessSet<micm::Matrix<double>, SparseMatrixTest>>(400, 100, 80);
}