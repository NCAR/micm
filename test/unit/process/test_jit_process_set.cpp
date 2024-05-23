#include "test_process_set_policy.hpp"

#include <micm/process/jit_process_set.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <random>

using Group2VectorMatrix = micm::VectorMatrix<double, 2>;
using Group200VectorMatrix = micm::VectorMatrix<double, 200>;
using Group300VectorMatrix = micm::VectorMatrix<double, 300>;
using Group400VectorMatrix = micm::VectorMatrix<double, 400>;

using Group2SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>;

TEST(JitProcessSet, VectorMatrix)
{
  testProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix, micm::JitProcessSet<2>>();
}

TEST(RandomJitProcessSet, VectorMatrix)
{
  testRandomSystem<Group200VectorMatrix, micm::StandardSparseMatrix, micm::JitProcessSet<200>>(200, 20, 30);
  testRandomSystem<Group300VectorMatrix, micm::StandardSparseMatrix, micm::JitProcessSet<300>>(300, 50, 40);
  testRandomSystem<Group300VectorMatrix, micm::StandardSparseMatrix, micm::JitProcessSet<300>>(300, 30, 20);
  testRandomSystem<Group400VectorMatrix, micm::StandardSparseMatrix, micm::JitProcessSet<400>>(400, 100, 80);
}