#include <gtest/gtest.h>

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include "test_sparse_matrix_policy.hpp"

TEST(SparseMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::SparseMatrixVectorOrdering<2>>();
}