#include <gtest/gtest.h>

#include <functional>
#include <micm/solver/linear_solver.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include "test_linear_solver_policy.hpp"

template<class T>
using SparseMatrixTest = micm::SparseMatrix<T>;

TEST(LinearSolver, DenseMatrixStandardOrdering)
{
  testDenseMatrix<micm::Matrix, SparseMatrixTest, micm::LinearSolver<double, SparseMatrixTest>>(
      [](const SparseMatrixTest<double>& matrix, double initial_value) -> micm::LinearSolver<double, SparseMatrixTest> {
        return micm::LinearSolver<double, SparseMatrixTest>{ matrix, initial_value };
      });
}

TEST(LinearSolver, RandomMatrixStandardOrdering)
{
  testRandomMatrix<micm::Matrix, SparseMatrixTest, micm::LinearSolver<double, SparseMatrixTest>>(
      [](const SparseMatrixTest<double>& matrix, double initial_value) -> micm::LinearSolver<double, SparseMatrixTest> {
        return micm::LinearSolver<double, SparseMatrixTest>{ matrix, initial_value };
      },
      5);
}

TEST(LinearSolver, DiagonalMatrixStandardOrdering)
{
  testDiagonalMatrix<micm::Matrix, SparseMatrixTest, micm::LinearSolver<double, SparseMatrixTest>>(
      [](const SparseMatrixTest<double>& matrix, double initial_value) -> micm::LinearSolver<double, SparseMatrixTest> {
        return micm::LinearSolver<double, SparseMatrixTest>{ matrix, initial_value };
      },
      5);
}

TEST(LinearSolver, DiagonalMarkowitzReorder)
{
  testMarkowitzReordering<micm::Matrix, SparseMatrixTest>();
}

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

TEST(LinearSolver, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<double, Group1SparseVectorMatrix>>(
      [](const Group1SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group1SparseVectorMatrix> {
        return micm::LinearSolver<double, Group1SparseVectorMatrix>{ matrix, initial_value };
      });
  testDenseMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<double, Group2SparseVectorMatrix>>(
      [](const Group2SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group2SparseVectorMatrix> {
        return micm::LinearSolver<double, Group2SparseVectorMatrix>{ matrix, initial_value };
      });
  testDenseMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<double, Group3SparseVectorMatrix>>(
      [](const Group3SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group3SparseVectorMatrix> {
        return micm::LinearSolver<double, Group3SparseVectorMatrix>{ matrix, initial_value };
      });
  testDenseMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<double, Group4SparseVectorMatrix>>(
      [](const Group4SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group4SparseVectorMatrix> {
        return micm::LinearSolver<double, Group4SparseVectorMatrix>{ matrix, initial_value };
      });
}

TEST(LinearSolver, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<double, Group1SparseVectorMatrix>>(
      [](const Group1SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group1SparseVectorMatrix> {
        return micm::LinearSolver<double, Group1SparseVectorMatrix>{ matrix, initial_value };
      },
      5);
  testRandomMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<double, Group2SparseVectorMatrix>>(
      [](const Group2SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group2SparseVectorMatrix> {
        return micm::LinearSolver<double, Group2SparseVectorMatrix>{ matrix, initial_value };
      },
      5);
  testRandomMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<double, Group3SparseVectorMatrix>>(
      [](const Group3SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group3SparseVectorMatrix> {
        return micm::LinearSolver<double, Group3SparseVectorMatrix>{ matrix, initial_value };
      },
      5);
  testRandomMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<double, Group4SparseVectorMatrix>>(
      [](const Group4SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group4SparseVectorMatrix> {
        return micm::LinearSolver<double, Group4SparseVectorMatrix>{ matrix, initial_value };
      },
      5);
}

TEST(LinearSolver, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::LinearSolver<double, Group1SparseVectorMatrix>>(
      [](const Group1SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group1SparseVectorMatrix> {
        return micm::LinearSolver<double, Group1SparseVectorMatrix>{ matrix, initial_value };
      },
      5);
  testDiagonalMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::LinearSolver<double, Group2SparseVectorMatrix>>(
      [](const Group2SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group2SparseVectorMatrix> {
        return micm::LinearSolver<double, Group2SparseVectorMatrix>{ matrix, initial_value };
      },
      5);
  testDiagonalMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::LinearSolver<double, Group3SparseVectorMatrix>>(
      [](const Group3SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group3SparseVectorMatrix> {
        return micm::LinearSolver<double, Group3SparseVectorMatrix>{ matrix, initial_value };
      },
      5);
  testDiagonalMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::LinearSolver<double, Group4SparseVectorMatrix>>(
      [](const Group4SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::LinearSolver<double, Group4SparseVectorMatrix> {
        return micm::LinearSolver<double, Group4SparseVectorMatrix>{ matrix, initial_value };
      },
      5);
}

TEST(LinearSolver, VectorDiagonalMarkowitzReordering)
{
  testMarkowitzReordering<Group1VectorMatrix, Group1SparseVectorMatrix>();
  testMarkowitzReordering<Group2VectorMatrix, Group2SparseVectorMatrix>();
  testMarkowitzReordering<Group3VectorMatrix, Group3SparseVectorMatrix>();
  testMarkowitzReordering<Group4VectorMatrix, Group4SparseVectorMatrix>();
}
