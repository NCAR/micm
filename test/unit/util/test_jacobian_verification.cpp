// Copyright (C) 2024-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/util/jacobian_verification.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering_compressed_sparse_row.hpp>

#include <gtest/gtest.h>

using DenseMatrix = micm::Matrix<double>;
using SparseMatrix = micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>;

// Simple 2-variable system: f(x,y) = [x*y, x^2 - y]
// Analytical Jacobian:
//   df0/dx = y,    df0/dy = x
//   df1/dx = 2*x,  df1/dy = -1
static void SimpleForcing(const DenseMatrix& vars, DenseMatrix& forcing)
{
  for (std::size_t block = 0; block < vars.NumRows(); ++block)
  {
    double x = vars[block][0];
    double y = vars[block][1];
    forcing[block][0] += x * y;
    forcing[block][1] += x * x - y;
  }
}

TEST(JacobianVerification, FiniteDifferenceMatchesAnalytical)
{
  const std::size_t num_species = 2;
  DenseMatrix variables(1, num_species, 0.0);
  variables[0][0] = 3.0;  // x
  variables[0][1] = 5.0;  // y

  auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(SimpleForcing, variables, num_species);

  // Expected: df0/dx = y = 5, df0/dy = x = 3, df1/dx = 2x = 6, df1/dy = -1
  EXPECT_NEAR(fd_jac[0][0 * num_species + 0], 5.0, 1e-5);   // df0/dx
  EXPECT_NEAR(fd_jac[0][0 * num_species + 1], 3.0, 1e-5);   // df0/dy
  EXPECT_NEAR(fd_jac[0][1 * num_species + 0], 6.0, 1e-5);   // df1/dx
  EXPECT_NEAR(fd_jac[0][1 * num_species + 1], -1.0, 1e-5);  // df1/dy
}

TEST(JacobianVerification, MultiBlockFiniteDifference)
{
  const std::size_t num_species = 2;
  DenseMatrix variables(3, num_species, 0.0);
  variables[0][0] = 1.0;
  variables[0][1] = 2.0;
  variables[1][0] = 0.5;
  variables[1][1] = 4.0;
  variables[2][0] = 10.0;
  variables[2][1] = 0.1;

  auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(SimpleForcing, variables, num_species);

  // Block 0: x=1, y=2
  EXPECT_NEAR(fd_jac[0][0], 2.0, 1e-5);   // df0/dx = y
  EXPECT_NEAR(fd_jac[0][1], 1.0, 1e-5);   // df0/dy = x
  EXPECT_NEAR(fd_jac[0][2], 2.0, 1e-5);   // df1/dx = 2x
  EXPECT_NEAR(fd_jac[0][3], -1.0, 1e-5);  // df1/dy

  // Block 1: x=0.5, y=4
  EXPECT_NEAR(fd_jac[1][0], 4.0, 1e-5);
  EXPECT_NEAR(fd_jac[1][1], 0.5, 1e-5);
  EXPECT_NEAR(fd_jac[1][2], 1.0, 1e-5);
  EXPECT_NEAR(fd_jac[1][3], -1.0, 1e-5);

  // Block 2: x=10, y=0.1
  EXPECT_NEAR(fd_jac[2][0], 0.1, 1e-5);
  EXPECT_NEAR(fd_jac[2][1], 10.0, 1e-5);
  EXPECT_NEAR(fd_jac[2][2], 20.0, 1e-5);
  EXPECT_NEAR(fd_jac[2][3], -1.0, 1e-5);
}

TEST(JacobianVerification, ComparePassesForCorrectJacobian)
{
  const std::size_t num_species = 2;
  DenseMatrix variables(1, num_species, 0.0);
  variables[0][0] = 3.0;
  variables[0][1] = 5.0;

  auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(SimpleForcing, variables, num_species);

  // Build sparse analytical Jacobian storing -df/dx (MICM convention)
  auto builder = SparseMatrix::Create(num_species)
                     .SetNumberOfBlocks(1)
                     .WithElement(0, 0)
                     .WithElement(0, 1)
                     .WithElement(1, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);
  SparseMatrix analytical{ builder };

  // Store -(df/dx): negate the true Jacobian values
  analytical[0][0][0] = -5.0;  // -df0/dx
  analytical[0][0][1] = -3.0;  // -df0/dy
  analytical[0][1][0] = -6.0;  // -df1/dx
  analytical[0][1][1] = 1.0;   // -df1/dy = -(-1) = 1

  auto result = micm::CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrix>(analytical, fd_jac, num_species);

  EXPECT_TRUE(result.passed);
  EXPECT_LT(result.max_abs_error, 1e-4);
}

TEST(JacobianVerification, CompareFailsForWrongJacobian)
{
  const std::size_t num_species = 2;
  DenseMatrix variables(1, num_species, 0.0);
  variables[0][0] = 3.0;
  variables[0][1] = 5.0;

  auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(SimpleForcing, variables, num_species);

  // Build an intentionally wrong Jacobian
  auto builder = SparseMatrix::Create(num_species)
                     .SetNumberOfBlocks(1)
                     .WithElement(0, 0)
                     .WithElement(0, 1)
                     .WithElement(1, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);
  SparseMatrix wrong_analytical{ builder };

  // Intentionally wrong values
  wrong_analytical[0][0][0] = -999.0;
  wrong_analytical[0][0][1] = -3.0;
  wrong_analytical[0][1][0] = -6.0;
  wrong_analytical[0][1][1] = 1.0;

  auto result = micm::CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrix>(wrong_analytical, fd_jac, num_species);

  EXPECT_FALSE(result.passed);
  EXPECT_EQ(result.worst_row, 0u);
  EXPECT_EQ(result.worst_col, 0u);
}

TEST(JacobianVerification, SparsityCompletenessPassesWhenComplete)
{
  const std::size_t num_species = 2;
  DenseMatrix variables(1, num_species, 0.0);
  variables[0][0] = 3.0;
  variables[0][1] = 5.0;

  auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(SimpleForcing, variables, num_species);

  // Declare all 4 entries as non-zero (fully dense)
  auto builder = SparseMatrix::Create(num_species)
                     .SetNumberOfBlocks(1)
                     .WithElement(0, 0)
                     .WithElement(0, 1)
                     .WithElement(1, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);
  SparseMatrix full_sparsity{ builder };

  auto result = micm::CheckJacobianSparsityCompleteness<DenseMatrix, SparseMatrix>(full_sparsity, fd_jac, num_species);

  EXPECT_TRUE(result.passed);
}

TEST(JacobianVerification, SparsityCompletenessFailsWhenMissingEntry)
{
  const std::size_t num_species = 2;
  DenseMatrix variables(1, num_species, 0.0);
  variables[0][0] = 3.0;
  variables[0][1] = 5.0;

  auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(SimpleForcing, variables, num_species);

  // Omit (1,0) from sparsity — but df1/dx = 2*x = 6 is non-zero
  auto builder = SparseMatrix::Create(num_species)
                     .SetNumberOfBlocks(1)
                     .WithElement(0, 0)
                     .WithElement(0, 1)
                     .WithElement(1, 1)
                     .InitialValue(0.0);
  SparseMatrix incomplete_sparsity{ builder };

  auto result = micm::CheckJacobianSparsityCompleteness<DenseMatrix, SparseMatrix>(incomplete_sparsity, fd_jac, num_species);

  EXPECT_FALSE(result.passed);
  EXPECT_EQ(result.worst_row, 1u);
  EXPECT_EQ(result.worst_col, 0u);
}

TEST(JacobianVerification, NearZeroVariableHandled)
{
  // Test that near-zero variables don't cause issues
  const std::size_t num_species = 2;
  DenseMatrix variables(1, num_species, 0.0);
  variables[0][0] = 0.0;  // exactly zero
  variables[0][1] = 1.0;

  auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(SimpleForcing, variables, num_species);

  // At x=0, y=1: df0/dx = y = 1, df0/dy = x = 0, df1/dx = 2x = 0, df1/dy = -1
  EXPECT_NEAR(fd_jac[0][0 * num_species + 0], 1.0, 1e-5);
  EXPECT_NEAR(fd_jac[0][0 * num_species + 1], 0.0, 1e-5);
  EXPECT_NEAR(fd_jac[0][1 * num_species + 0], 0.0, 1e-5);
  EXPECT_NEAR(fd_jac[0][1 * num_species + 1], -1.0, 1e-5);
}
