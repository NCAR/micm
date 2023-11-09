#pragma once
#include <gtest/gtest.h>
#include <functional>
#include <random>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/solver/cuda_linear_solver.hpp>
#include <micm/solver/cuda_lu_decomposition.hpp>
#include "test_linear_solver_policy.hpp"

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

TEST(CudaLinearSolver, DenseMatrixVectorOrdering)
{
  testDenseMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::CudaLinearSolver<double, Group1SparseVectorMatrix, micm::CudaLuDecomposition>>(
      [](const Group1SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group1SparseVectorMatrix, micm::CudaLuDecomposition> {
        return micm::CudaLinearSolver<double, Group1SparseVectorMatrix, micm::CudaLuDecomposition>{matrix, initial_value};
      });
}

TEST(CudaLinearSolver, RandomMatrixVectorOrdering)
{
  testRandomMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::CudaLinearSolver<double, Group1SparseVectorMatrix>>(
      [](const Group1SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group1SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group1SparseVectorMatrix>{ matrix, initial_value };
      },
      1);
  testRandomMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::CudaLinearSolver<double, Group2SparseVectorMatrix>>(
      [](const Group2SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group2SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group2SparseVectorMatrix>{ matrix, initial_value };
      },
      2);
  testRandomMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::CudaLinearSolver<double, Group3SparseVectorMatrix>>(
      [](const Group3SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group3SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group3SparseVectorMatrix>{ matrix, initial_value };
      },
      3);
  testRandomMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::CudaLinearSolver<double, Group4SparseVectorMatrix>>(
      [](const Group4SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group4SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group4SparseVectorMatrix>{ matrix, initial_value };
      },
      4);
}

TEST(CudaLinearSolver, DiagonalMatrixVectorOrdering)
{
  testDiagonalMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::CudaLinearSolver<double, Group1SparseVectorMatrix>>(
      [](const Group1SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group1SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group1SparseVectorMatrix>{ matrix, initial_value };
      },
      1);
  testDiagonalMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::CudaLinearSolver<double, Group2SparseVectorMatrix>>(
      [](const Group2SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group2SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group2SparseVectorMatrix>{ matrix, initial_value };
      },
      2);
  testDiagonalMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::CudaLinearSolver<double, Group3SparseVectorMatrix>>(
      [](const Group3SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group3SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group3SparseVectorMatrix>{ matrix, initial_value };
      },
      3);
  testDiagonalMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::CudaLinearSolver<double, Group4SparseVectorMatrix>>(
      [](const Group4SparseVectorMatrix<double>& matrix,
         double initial_value) -> micm::CudaLinearSolver<double, Group4SparseVectorMatrix> {
        return micm::CudaLinearSolver<double, Group4SparseVectorMatrix>{ matrix, initial_value };
      },
      4);
}


