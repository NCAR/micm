#include "../../solver/test_lu_decomposition_in_place_policy.hpp"

#include <micm/cuda/solver/cuda_lu_decomposition_mozart_in_place.hpp>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/solver/lu_decomposition_mozart_in_place.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <random>
#include <vector>

using Group1CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group100CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<100>>;
using Group1000CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1000>>;
using Group100000CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<100000>>;

TEST(CudaLuDecompositionMozartInPlace, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(1);
  testRandomMatrix<Group100CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(100);
  testRandomMatrix<Group1000CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(1000);
  testRandomMatrix<Group100000CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(100000);
}

TEST(CudaLuDecompositionMozartInPlace, AgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    testExtremeValueInitialization<Group1CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(1, value);
    testExtremeValueInitialization<Group100CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(100, value);
    testExtremeValueInitialization<Group1000CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(1000, value);
    testExtremeValueInitialization<Group100000CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(100000, value);
  }
}