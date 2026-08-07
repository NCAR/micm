#include "test_lu_decomposition_in_place_policy.hpp"

#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

using SparseMatrixTest = micm::SparseMatrix<micm::Real>;

using Group1SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<4>>;

using LUStandard = micm::LuDecompositionDoolittleInPlace<SparseMatrixTest>;
using LUVector1  = micm::LuDecompositionDoolittleInPlace<Group1SparseVectorMatrix>;
using LUVector2  = micm::LuDecompositionDoolittleInPlace<Group2SparseVectorMatrix>;
using LUVector3  = micm::LuDecompositionDoolittleInPlace<Group3SparseVectorMatrix>;
using LUVector4  = micm::LuDecompositionDoolittleInPlace<Group4SparseVectorMatrix>;

TEST(LuDecompositionDoolittleInPlace, DenseMatrixStandardOrdering)
{
  TestDenseMatrix<SparseMatrixTest, LUStandard>();
}

TEST(LuDecompositionDoolittleInPlace, RandomMatrixStandardOrdering)
{
  TestRandomMatrix<SparseMatrixTest, LUStandard>(1);
  TestRandomMatrix<SparseMatrixTest, LUStandard>(5);
}

TEST(LuDecompositionDoolittleInPlace, DiagonalMatrixStandardOrdering)
{
  TestDiagonalMatrix<SparseMatrixTest, LUStandard>(5);
}

TEST(LuDecompositionDoolittleInPlace, AgnosticToInitialValueStandardOrdering)
{
  micm::Real initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<SparseMatrixTest, LUStandard>(5, value);
  }
}

TEST(LuDecompositionDoolittleInPlace, DenseMatrixVectorOrdering)
{
  TestDenseMatrix<Group1SparseVectorMatrix, LUVector1>();
  TestDenseMatrix<Group2SparseVectorMatrix, LUVector2>();
  TestDenseMatrix<Group3SparseVectorMatrix, LUVector3>();
  TestDenseMatrix<Group4SparseVectorMatrix, LUVector4>();
}

TEST(LuDecompositionDoolittleInPlace, RandomMatrixVectorOrdering)
{
  TestRandomMatrix<Group1SparseVectorMatrix, LUVector1>(5);
  TestRandomMatrix<Group2SparseVectorMatrix, LUVector2>(5);
  TestRandomMatrix<Group3SparseVectorMatrix, LUVector3>(5);
  TestRandomMatrix<Group4SparseVectorMatrix, LUVector4>(5);
}

TEST(LuDecompositionDoolittleInPlace, DiagonalMatrixVectorOrdering)
{
  TestDiagonalMatrix<Group1SparseVectorMatrix, LUVector1>(5);
  TestDiagonalMatrix<Group2SparseVectorMatrix, LUVector2>(5);
  TestDiagonalMatrix<Group3SparseVectorMatrix, LUVector3>(5);
  TestDiagonalMatrix<Group4SparseVectorMatrix, LUVector4>(5);
}

TEST(LuDecompositionDoolittleInPlace, AgnosticToInitialValueVectorOrdering)
{
  micm::Real initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<Group1SparseVectorMatrix, LUVector1>(5, value);
    TestExtremeValueInitialization<Group2SparseVectorMatrix, LUVector2>(5, value);
    TestExtremeValueInitialization<Group3SparseVectorMatrix, LUVector3>(5, value);
    TestExtremeValueInitialization<Group4SparseVectorMatrix, LUVector4>(5, value);
  }
}