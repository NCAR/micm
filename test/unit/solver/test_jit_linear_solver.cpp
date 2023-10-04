#include <gtest/gtest.h>

#include <functional>
#include <micm/solver/jit_linear_solver.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

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

TEST(JitLinearSolver, DenseMatrixVectorOrdering)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  testDenseMatrix<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::JitLinearSolver<1, Group1SparseVectorMatrix, micm::JitLuDecomposition<1>>>(
      [&](const Group1SparseVectorMatrix<double>& matrix,
          double initial_value) -> micm::JitLinearSolver<1, Group1SparseVectorMatrix, micm::JitLuDecomposition<1>>
      {
        return micm::JitLinearSolver<1, Group1SparseVectorMatrix, micm::JitLuDecomposition<1>>{ jit.get(),
                                                                                                matrix,
                                                                                                initial_value };
      });
}

TEST(JitLinearSolver, RandomMatrixVectorOrdering)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  testRandomMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::JitLinearSolver<1, Group1SparseVectorMatrix>>(
      [&](const Group1SparseVectorMatrix<double>& matrix,
          double initial_value) -> micm::JitLinearSolver<1, Group1SparseVectorMatrix> {
        return micm::JitLinearSolver<1, Group1SparseVectorMatrix>{ jit.get(), matrix, initial_value };
      },
      1);
  testRandomMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::JitLinearSolver<2, Group2SparseVectorMatrix>>(
      [&](const Group2SparseVectorMatrix<double>& matrix,
          double initial_value) -> micm::JitLinearSolver<2, Group2SparseVectorMatrix> {
        return micm::JitLinearSolver<2, Group2SparseVectorMatrix>{ jit.get(), matrix, initial_value };
      },
      2);
  testRandomMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::JitLinearSolver<3, Group3SparseVectorMatrix>>(
      [&](const Group3SparseVectorMatrix<double>& matrix,
          double initial_value) -> micm::JitLinearSolver<3, Group3SparseVectorMatrix> {
        return micm::JitLinearSolver<3, Group3SparseVectorMatrix>{ jit.get(), matrix, initial_value };
      },
      3);
  testRandomMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::JitLinearSolver<4, Group4SparseVectorMatrix>>(
      [&](const Group4SparseVectorMatrix<double>& matrix,
          double initial_value) -> micm::JitLinearSolver<4, Group4SparseVectorMatrix> {
        return micm::JitLinearSolver<4, Group4SparseVectorMatrix>{ jit.get(), matrix, initial_value };
      },
      4);
}

TEST(JitLinearSolver, DiagonalMatrixVectorOrdering)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  testDiagonalMatrix<Group1VectorMatrix, Group1SparseVectorMatrix, micm::JitLinearSolver<1, Group1SparseVectorMatrix>>(
      [&](const Group1SparseVectorMatrix<double>& matrix,
          double initial_value) -> micm::JitLinearSolver<1, Group1SparseVectorMatrix> {
        return micm::JitLinearSolver<1, Group1SparseVectorMatrix>{ jit.get(), matrix, initial_value };
      },
      1);
  testDiagonalMatrix<Group2VectorMatrix, Group2SparseVectorMatrix, micm::JitLinearSolver<2, Group2SparseVectorMatrix>>(
      [&](const Group2SparseVectorMatrix<double>& matrix,
          double initial_value) -> micm::JitLinearSolver<2, Group2SparseVectorMatrix> {
        return micm::JitLinearSolver<2, Group2SparseVectorMatrix>{ jit.get(), matrix, initial_value };
      },
      2);
  testDiagonalMatrix<Group3VectorMatrix, Group3SparseVectorMatrix, micm::JitLinearSolver<3, Group3SparseVectorMatrix>>(
      [&](const Group3SparseVectorMatrix<double>& matrix,
          double initial_value) -> micm::JitLinearSolver<3, Group3SparseVectorMatrix> {
        return micm::JitLinearSolver<3, Group3SparseVectorMatrix>{ jit.get(), matrix, initial_value };
      },
      3);
  testDiagonalMatrix<Group4VectorMatrix, Group4SparseVectorMatrix, micm::JitLinearSolver<4, Group4SparseVectorMatrix>>(
      [&](const Group4SparseVectorMatrix<double>& matrix,
          double initial_value) -> micm::JitLinearSolver<4, Group4SparseVectorMatrix> {
        return micm::JitLinearSolver<4, Group4SparseVectorMatrix>{ jit.get(), matrix, initial_value };
      },
      4);
}