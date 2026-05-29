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
using Group3CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group27CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<27>>;
using Group32CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<32>>;
using Group43CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<43>>;
using Group77CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<77>>;
using Group113CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<113>>;
using Group193CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<193>>;
using Group281CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<281>>;
using Group472CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<472>>;
using Group512CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<512>>;
using Group739CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<739>>;
using Group1130CudaSparseMatrix = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1130>>;

TEST(CudaLuDecompositionMozartInPlace, RandomMatrixVectorOrdering)
{
  TestRandomMatrix<Group1CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group3CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group27CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group32CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group43CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group77CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group113CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group193CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group281CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group472CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group512CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group739CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
  TestRandomMatrix<Group1130CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400);
}

TEST(CudaLuDecompositionMozartInPlace, AgnosticToInitialValue)
{
  double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<Group1CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group3CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group27CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group32CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group43CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group77CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group113CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group193CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group281CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group472CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group512CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group739CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
    TestExtremeValueInitialization<Group1130CudaSparseMatrix, micm::CudaLuDecompositionMozartInPlace>(400, value);
  }
}