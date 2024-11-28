#include "../../solver/test_lu_decomposition_policy.hpp"

#include <micm/jit/solver/jit_lu_decomposition_doolittle.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

using Group1SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<4>>;

TEST(JitLuDecompositionDoolittle, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1SparseVectorMatrix, micm::JitLuDecompositionDoolittle<1>>();
  testDenseMatrix<Group2SparseVectorMatrix, micm::JitLuDecompositionDoolittle<2>>();
  testDenseMatrix<Group3SparseVectorMatrix, micm::JitLuDecompositionDoolittle<3>>();
  testDenseMatrix<Group4SparseVectorMatrix, micm::JitLuDecompositionDoolittle<4>>();
}

TEST(JitLuDecompositionDoolittle, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1SparseVectorMatrix, micm::JitLuDecompositionDoolittle<1>>(1);
  testRandomMatrix<Group2SparseVectorMatrix, micm::JitLuDecompositionDoolittle<2>>(2);
  testRandomMatrix<Group3SparseVectorMatrix, micm::JitLuDecompositionDoolittle<3>>(3);
  testRandomMatrix<Group4SparseVectorMatrix, micm::JitLuDecompositionDoolittle<4>>(4);
}

TEST(JitLuDecompositionDoolittle, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1SparseVectorMatrix, micm::JitLuDecompositionDoolittle<1>>(1);
  testDiagonalMatrix<Group2SparseVectorMatrix, micm::JitLuDecompositionDoolittle<2>>(2);
  testDiagonalMatrix<Group3SparseVectorMatrix, micm::JitLuDecompositionDoolittle<3>>(3);
  testDiagonalMatrix<Group4SparseVectorMatrix, micm::JitLuDecompositionDoolittle<4>>(4);
}

TEST(JitLuDecompositionDoolittle, AgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<Group1SparseVectorMatrix, micm::JitLuDecompositionDoolittle<1>>(1, value);
    testExtremeValueInitialization<Group2SparseVectorMatrix, micm::JitLuDecompositionDoolittle<2>>(2, value);
    testExtremeValueInitialization<Group3SparseVectorMatrix, micm::JitLuDecompositionDoolittle<3>>(3, value);
    testExtremeValueInitialization<Group4SparseVectorMatrix, micm::JitLuDecompositionDoolittle<4>>(4, value);
  }
}