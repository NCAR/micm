#include <gtest/gtest.h>

#include <micm/solver/jit_lu_decomposition.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include "test_lu_decomposition_policy.hpp"

template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;

TEST(JitLuDecomposition, DenseMatrixVectorOrdering)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  testDenseMatrix<Group1SparseVectorMatrix, micm::JitLuDecomposition<1>>(
      [&](const Group1SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<1> {
        return micm::JitLuDecomposition<1>{ jit.get(), matrix };
      });
  testDenseMatrix<Group2SparseVectorMatrix, micm::JitLuDecomposition<2>>(
      [&](const Group2SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<2> {
        return micm::JitLuDecomposition<2>{ jit.get(), matrix };
      });
  testDenseMatrix<Group3SparseVectorMatrix, micm::JitLuDecomposition<3>>(
      [&](const Group3SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<3> {
        return micm::JitLuDecomposition<3>{ jit.get(), matrix };
      });
  testDenseMatrix<Group4SparseVectorMatrix, micm::JitLuDecomposition<4>>(
      [&](const Group4SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<4> {
        return micm::JitLuDecomposition<4>{ jit.get(), matrix };
      });
}

TEST(JitLuDecomposition, SingularMatrixVectorOrdering)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  testSingularMatrix<Group1SparseVectorMatrix, micm::JitLuDecomposition<1>>(
      [&](const Group1SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<1> {
        return micm::JitLuDecomposition<1>{ jit.get(), matrix };
      });
  testSingularMatrix<Group2SparseVectorMatrix, micm::JitLuDecomposition<2>>(
      [&](const Group2SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<2> {
        return micm::JitLuDecomposition<2>{ jit.get(), matrix };
      });
  testSingularMatrix<Group3SparseVectorMatrix, micm::JitLuDecomposition<3>>(
      [&](const Group3SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<3> {
        return micm::JitLuDecomposition<3>{ jit.get(), matrix };
      });
  testSingularMatrix<Group4SparseVectorMatrix, micm::JitLuDecomposition<4>>(
      [&](const Group4SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<4> {
        return micm::JitLuDecomposition<4>{ jit.get(), matrix };
      });
}

TEST(JitLuDecomposition, RandomMatrixVectorOrdering)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  testRandomMatrix<Group1SparseVectorMatrix, micm::JitLuDecomposition<1>>(
      [&](const Group1SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<1> {
        return micm::JitLuDecomposition<1>{ jit.get(), matrix };
      },
      1);
  testRandomMatrix<Group2SparseVectorMatrix, micm::JitLuDecomposition<2>>(
      [&](const Group2SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<2> {
        return micm::JitLuDecomposition<2>{ jit.get(), matrix };
      },
      2);
  testRandomMatrix<Group3SparseVectorMatrix, micm::JitLuDecomposition<3>>(
      [&](const Group3SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<3> {
        return micm::JitLuDecomposition<3>{ jit.get(), matrix };
      },
      3);
  testRandomMatrix<Group4SparseVectorMatrix, micm::JitLuDecomposition<4>>(
      [&](const Group4SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<4> {
        return micm::JitLuDecomposition<4>{ jit.get(), matrix };
      },
      4);
}

TEST(JitLuDecomposition, DiagonalMatrixVectorOrdering)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  testDiagonalMatrix<Group1SparseVectorMatrix, micm::JitLuDecomposition<1>>(
      [&](const Group1SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<1> {
        return micm::JitLuDecomposition<1>{ jit.get(), matrix };
      },
      1);
  testDiagonalMatrix<Group2SparseVectorMatrix, micm::JitLuDecomposition<2>>(
      [&](const Group2SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<2> {
        return micm::JitLuDecomposition<2>{ jit.get(), matrix };
      },
      2);
  testDiagonalMatrix<Group3SparseVectorMatrix, micm::JitLuDecomposition<3>>(
      [&](const Group3SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<3> {
        return micm::JitLuDecomposition<3>{ jit.get(), matrix };
      },
      3);
  testDiagonalMatrix<Group4SparseVectorMatrix, micm::JitLuDecomposition<4>>(
      [&](const Group4SparseVectorMatrix<double>& matrix) -> micm::JitLuDecomposition<4> {
        return micm::JitLuDecomposition<4>{ jit.get(), matrix };
      },
      4);
}
