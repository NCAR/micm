#include <gtest/gtest.h>

#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include "test_lu_decomposition_policy.hpp"

template<class T>
using SparseMatrixTest = micm::SparseMatrix<T>;

template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;

TEST(LuDecomposition, DenseMatrixStandardOrdering)
{
  testDenseMatrix<SparseMatrixTest, micm::LuDecomposition>(
      [](const SparseMatrixTest<double>& matrix) -> micm::LuDecomposition { return micm::LuDecomposition{ matrix }; });
}

TEST(LuDecomposition, SingularMatrixStandardOrdering)
{
  testSingularMatrix<SparseMatrixTest, micm::LuDecomposition>(
      [](const SparseMatrixTest<double>& matrix) -> micm::LuDecomposition { return micm::LuDecomposition{ matrix }; });
}

TEST(LuDecomposition, RandomMatrixStandardOrdering)
{
  testRandomMatrix<SparseMatrixTest, micm::LuDecomposition>(
      [](const SparseMatrixTest<double>& matrix) -> micm::LuDecomposition { return micm::LuDecomposition{ matrix }; }, 5);
}

TEST(LuDecomposition, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<SparseMatrixTest, micm::LuDecomposition>(
      [](const SparseMatrixTest<double>& matrix) -> micm::LuDecomposition { return micm::LuDecomposition{ matrix }; }, 5);
}

TEST(LuDecomposition, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group1SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; });
  testDenseMatrix<Group2SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group2SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; });
  testDenseMatrix<Group3SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group3SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; });
  testDenseMatrix<Group4SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group4SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; });
}

TEST(LuDecomposition, SingluarMatrixVectorOrdering)
{
  testSingularMatrix<Group1SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group1SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; });
  testSingularMatrix<Group2SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group2SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; });
  testSingularMatrix<Group3SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group3SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; });
  testSingularMatrix<Group4SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group4SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; });
}

TEST(LuDecomposition, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group1SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; },
      5);
  testRandomMatrix<Group2SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group2SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; },
      5);
  testRandomMatrix<Group3SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group3SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; },
      5);
  testRandomMatrix<Group4SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group4SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; },
      5);
}

TEST(LuDecomposition, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group1SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; },
      5);
  testDiagonalMatrix<Group2SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group2SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; },
      5);
  testDiagonalMatrix<Group3SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group3SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; },
      5);
  testDiagonalMatrix<Group4SparseVectorMatrix, micm::LuDecomposition>(
      [](const Group4SparseVectorMatrix<double>& matrix) -> micm::LuDecomposition
      { return micm::LuDecomposition{ matrix }; },
      5);
}
