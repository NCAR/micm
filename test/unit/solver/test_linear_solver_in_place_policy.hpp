#include <micm/solver/linear_solver.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <random>
#include <type_traits>

// The comparator f receives the expected value, the reconstructed A*x value, and a scale equal to the
// sum of the magnitudes of the terms that were summed to form that row's result. The scale lets callers
// express a tolerance relative to the size of the arithmetic (i.e. the backward error of the solve),
// which is the meaningful bound when the result is formed by cancellation of much larger terms.
template<typename T, class MatrixPolicy, class SparseMatrixPolicy>
void CheckResults(
    const SparseMatrixPolicy& A,
    const MatrixPolicy& b,
    const MatrixPolicy& x,
    const std::function<void(const T, const T, const T)>& f)
{
  T result;
  EXPECT_EQ(A.NumberOfBlocks(), b.NumRows());
  EXPECT_EQ(A.NumberOfBlocks(), x.NumRows());
  for (micm::Index i_block = 0; i_block < A.NumberOfBlocks(); ++i_block)
  {
    for (micm::Index i = 0; i < A.NumRows(); ++i)
    {
      result = 0.0;
      T scale = 0.0;
      for (micm::Index j = 0; j < A.NumColumns(); ++j)
      {
        if (!A.IsZero(i, j))
        {
          result += A[i_block][i][j] * x[i_block][j];
          scale += std::abs(A[i_block][i][j] * x[i_block][j]);
        }
      }
      f(b[i_block][i], result, scale);
    }
  }
}

template<class SparseMatrixPolicy>
void PrintMatrix(const SparseMatrixPolicy& matrix, micm::Index width)
{
  for (micm::Index i_block = 0; i_block < matrix.NumberOfBlocks(); ++i_block)
  {
    std::cout << "block: " << i_block << std::endl;
    for (micm::Index i = 0; i < matrix.NumRows(); ++i)
    {
      for (micm::Index j = 0; j < matrix.NumColumns(); ++j)
      {
        if (matrix.IsZero(i, j))
        {
          std::cout << " " << std::setfill('-') << std::setw(width) << "-" << " ";
        }
        else
        {
          std::cout << " " << std::setfill(' ') << std::setw(width) << matrix[i_block][i][j] << " ";
        }
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

template<class MatrixPolicy, class SparseMatrixPolicy, class LinearSolverPolicy>
void TestDenseMatrix()
{
  using FloatingPointType = typename MatrixPolicy::value_type;

  SparseMatrixPolicy A = SparseMatrixPolicy(SparseMatrixPolicy::Create(3)
                                                .InitialValue(0)
                                                .WithElement(0, 0)
                                                .WithElement(0, 1)
                                                .WithElement(0, 2)
                                                .WithElement(1, 0)
                                                .WithElement(1, 1)
                                                .WithElement(1, 2)
                                                .WithElement(2, 0)
                                                .WithElement(2, 1)
                                                .WithElement(2, 2));
  MatrixPolicy b(1, 3, 0.0);
  MatrixPolicy x(1, 3, 0.0);

  A[0][0][0] = 2;
  A[0][0][1] = -1;
  A[0][0][2] = -2;
  A[0][1][0] = -4;
  A[0][1][1] = 6;
  A[0][1][2] = 3;
  A[0][2][0] = -4;
  A[0][2][1] = -2;
  A[0][2][2] = 8;

  b[0][0] = 23;
  b[0][1] = 42;
  b[0][2] = 9;

  x = b;

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<MatrixPolicy>(b);
  CheckCopyToDevice<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(A, 0);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(A, 0, false);
  alu.Fill(0);
  for (micm::Index i = 0; i < A.NumRows(); ++i)
  {
    for (micm::Index j = 0; j < A.NumColumns(); ++j)
    {
      if (!A.IsZero(i, j))
      {
        alu[0][i][j] = A[0][i][j];
      }
    }
  }

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<SparseMatrixPolicy>(alu);

  solver.Factor(alu);
  solver.template Solve<MatrixPolicy>(x, alu);

  // Only copy the data to the host when it is a CudaMatrix
  CheckCopyToHost<MatrixPolicy>(x);

  // In FLOAT precision the reconstructed b is only accurate to the backward error of the solve, which scales
  // with the magnitude of the summed terms; keep DOUBLE at its original absolute bound.
  CheckResults<FloatingPointType, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const FloatingPointType a, const FloatingPointType b, const FloatingPointType scale) -> void {
        const FloatingPointType tol =
            std::is_same_v<micm::Real, double> ? 1.0e-5 : scale * 1.0e5 * std::numeric_limits<micm::Real>::epsilon();
        EXPECT_NEAR(a, b, tol);
      });
}

template<class MatrixPolicy, class SparseMatrixPolicy, class LinearSolverPolicy>
void TestRandomMatrix(micm::Index number_of_blocks)
{
  using FloatingPointType = typename MatrixPolicy::value_type;

  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy::Create(10).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (micm::Index i = 0; i < 10; ++i)
  {
    for (micm::Index j = 0; j < 10; ++j)
    {
      if (i == j || gen_bool())
      {
        builder = builder.WithElement(i, j);
      }
    }
  }

  SparseMatrixPolicy A(builder);
  MatrixPolicy b(number_of_blocks, 10, 0.0);
  MatrixPolicy x(number_of_blocks, 10, 0.0);

  for (micm::Index i = 0; i < 10; ++i)
  {
    for (micm::Index j = 0; j < 10; ++j)
    {
      if (!A.IsZero(i, j))
      {
        for (micm::Index i_block = 0; i_block < number_of_blocks; ++i_block)
        {
          A[i_block][i][j] = get_double();
        }
      }
    }
  }

  for (micm::Index i = 0; i < 10; ++i)
  {
    for (micm::Index i_block = 0; i_block < number_of_blocks; ++i_block)
    {
      b[i_block][i] = get_double();
    }
  }

  x = b;

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(A, 0);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(A, 0, false);
  alu.Fill(0);
  for (micm::Index i_block = 0; i_block < number_of_blocks; ++i_block)
  {
    for (micm::Index i = 0; i < A.NumRows(); ++i)
    {
      for (micm::Index j = 0; j < A.NumColumns(); ++j)
      {
        if (!A.IsZero(i, j))
        {
          alu[i_block][i][j] = A[i_block][i][j];
        }
      }
    }
  }

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<SparseMatrixPolicy>(alu);

  solver.Factor(alu);
  solver.template Solve<MatrixPolicy>(x, alu);

  // Only copy the data to the host when it is a CudaMatrix
  CheckCopyToHost<MatrixPolicy>(x);

  // In FLOAT precision the reconstructed b is only accurate to the backward error of the solve, which scales
  // with the magnitude of the summed terms; keep DOUBLE at its original absolute bound.
  CheckResults<FloatingPointType, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const FloatingPointType a, const FloatingPointType b, const FloatingPointType scale) -> void {
        const FloatingPointType tol =
            std::is_same_v<micm::Real, double> ? 1.0e-6 : scale * 1.0e5 * std::numeric_limits<micm::Real>::epsilon();
        EXPECT_NEAR(a, b, tol);
      });
}

template<class MatrixPolicy, class SparseMatrixPolicy, class LinearSolverPolicy>
void TestExtremeInitialValue(micm::Index number_of_blocks, micm::Real initial_value)
{
  using FloatingPointType = typename MatrixPolicy::value_type;

  const unsigned int seed = 12345;
  std::mt19937 generator(seed);
  const micm::Real point_five = 0.5;
  const micm::Real two = 2.0;

  auto gen_bool = std::bind(std::bernoulli_distribution(point_five), generator);
  auto get_double = std::bind(std::lognormal_distribution<micm::Real>(-two, two), generator);
  const micm::Index size = 30;

  auto builder = SparseMatrixPolicy::Create(size).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (micm::Index i = 0; i < size; ++i)
  {
    for (micm::Index j = 0; j < size; ++j)
    {
      if (i == j || gen_bool())
      {
        builder = builder.WithElement(i, j);
      }
    }
  }

  SparseMatrixPolicy A(builder);
  MatrixPolicy b(number_of_blocks, size, 0.0);
  MatrixPolicy x(number_of_blocks, size, 0.0);

  for (micm::Index i = 0; i < size; ++i)
  {
    for (micm::Index j = 0; j < size; ++j)
    {
      if (!A.IsZero(i, j))
      {
        for (micm::Index i_block = 0; i_block < number_of_blocks; ++i_block)
        {
          A[i_block][i][j] = get_double();
        }
      }
    }
  }

  for (micm::Index i = 0; i < size; ++i)
  {
    for (micm::Index i_block = 0; i_block < number_of_blocks; ++i_block)
    {
      b[i_block][i] = get_double();
    }
  }

  x = b;

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(A, initial_value);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(A, initial_value, false);
  for (micm::Index i_block = 0; i_block < number_of_blocks; ++i_block)
  {
    for (micm::Index i = 0; i < A.NumRows(); ++i)
    {
      for (micm::Index j = 0; j < A.NumColumns(); ++j)
      {
        if (!A.IsZero(i, j))
        {
          alu[i_block][i][j] = A[i_block][i][j];
        }
        else if (!alu.IsZero(i, j))
        {
          alu[i_block][i][j] = 0;
        }
      }
    }
  }

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<SparseMatrixPolicy>(alu);

  solver.Factor(alu);

  // Only copy the data to the host when it is a CudaMatrix
  CheckCopyToHost<SparseMatrixPolicy>(alu);

  solver.template Solve<MatrixPolicy>(x, alu);

  // Only copy the data to the host when it is a CudaMatrix
  CheckCopyToHost<MatrixPolicy>(x);

  // In FLOAT precision the reconstructed b is only accurate to the backward error of the solve, which scales
  // with the magnitude of the summed terms; keep DOUBLE at its original absolute bound.
  CheckResults<FloatingPointType, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const FloatingPointType a, const FloatingPointType b, const FloatingPointType scale) -> void {
        const FloatingPointType tol =
            std::is_same_v<micm::Real, double> ? 2.0e-06 : scale * 1.0e5 * std::numeric_limits<micm::Real>::epsilon();
        EXPECT_NEAR(a, b, tol);
      });
}

template<class MatrixPolicy, class SparseMatrixPolicy, class LinearSolverPolicy>
void TestDiagonalMatrix(micm::Index number_of_blocks)
{
  using FloatingPointType = typename MatrixPolicy::value_type;

  auto get_double = std::bind(std::lognormal_distribution(-2.0, 4.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy::Create(6).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (micm::Index i = 0; i < 6; ++i)
  {
    builder = builder.WithElement(i, i);
  }

  SparseMatrixPolicy A(builder);
  MatrixPolicy b(number_of_blocks, 6, 0.0);
  MatrixPolicy x(number_of_blocks, 6, 0.0);

  for (micm::Index i = 0; i < 6; ++i)
  {
    for (micm::Index i_block = 0; i_block < number_of_blocks; ++i_block)
    {
      A[i_block][i][i] = get_double();
    }
  }

  for (micm::Index i = 0; i < 6; ++i)
  {
    for (micm::Index i_block = 0; i_block < number_of_blocks; ++i_block)
    {
      b[i_block][i] = get_double();
    }
  }

  x = b;

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(A, 0);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(A, 0, false);
  alu.Fill(0);
  for (micm::Index i_block = 0; i_block < number_of_blocks; ++i_block)
  {
    for (micm::Index i = 0; i < A.NumRows(); ++i)
    {
      for (micm::Index j = 0; j < A.NumColumns(); ++j)
      {
        if (!A.IsZero(i, j))
        {
          alu[i_block][i][j] = A[i_block][i][j];
        }
      }
    }
  }

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<SparseMatrixPolicy>(alu);

  solver.Factor(alu);
  solver.template Solve<MatrixPolicy>(x, alu);

  // Only copy the data to the host when it is a CudaMatrix
  CheckCopyToHost<MatrixPolicy>(x);

  CheckResults<FloatingPointType, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const FloatingPointType a, const FloatingPointType b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}

template<class MatrixPolicy, class SparseMatrixPolicy>
void TestMarkowitzReordering()
{
  const micm::Index order = 50;
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  MatrixPolicy orig(order, order, 0);

  for (micm::Index i = 0; i < order; ++i)
  {
    for (micm::Index j = 0; j < order; ++j)
    {
      orig[i][j] = (i == j || gen_bool()) ? 1 : 0;
    }
  }

  auto reorder_map = micm::DiagonalMarkowitzReorder<MatrixPolicy>(orig);

  auto builder = SparseMatrixPolicy::Create(50);
  for (micm::Index i = 0; i < order; ++i)
  {
    for (micm::Index j = 0; j < order; ++j)
    {
      if (orig[i][j] != 0)
      {
        builder = builder.WithElement(i, j);
      }
    }
  }
  SparseMatrixPolicy orig_jac{ builder };

  builder = SparseMatrixPolicy::Create(50);
  for (micm::Index i = 0; i < order; ++i)
  {
    for (micm::Index j = 0; j < order; ++j)
    {
      if (orig[reorder_map[i]][reorder_map[j]] != 0)
      {
        builder = builder.WithElement(i, j);
      }
    }
  }
  SparseMatrixPolicy reordered_jac{ builder };

  auto orig_LU_calc = micm::LuDecompositionInPlace::Create<SparseMatrixPolicy>(orig_jac);
  auto reordered_LU_calc = micm::LuDecompositionInPlace::Create<SparseMatrixPolicy>(reordered_jac);

  auto orig_LU = orig_LU_calc.template GetLUMatrix<SparseMatrixPolicy>(orig_jac, 0.0, false);
  auto reordered_LU = reordered_LU_calc.template GetLUMatrix<SparseMatrixPolicy>(reordered_jac, 0.0, false);

  micm::Index sum_orig = 0;
  micm::Index sum_reordered = 0;
  for (micm::Index i = 0; i < reorder_map.size(); ++i)
  {
    sum_orig += i;
    sum_reordered += reorder_map[i];
  }

  EXPECT_EQ(sum_orig, sum_reordered);
  EXPECT_GT(orig_LU.FlatBlockSize(), reordered_LU.FlatBlockSize());
}
