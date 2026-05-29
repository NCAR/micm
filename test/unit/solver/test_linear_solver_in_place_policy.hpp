#include <micm/solver/linear_solver.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

#include <functional>
#include <iomanip>
#include <random>

template<typename T, class MatrixPolicy, class SparseMatrixPolicy>
void CheckResults(
    const SparseMatrixPolicy A,
    const MatrixPolicy B,
    const MatrixPolicy X,
    const std::function<void(const T, const T)> F)
{
  T result;
  EXPECT_EQ(A.NumberOfBlocks(), B.NumRows());
  EXPECT_EQ(A.NumberOfBlocks(), X.NumRows());
  for (std::size_t i_block = 0; i_block < A.NumberOfBlocks(); ++i_block)
  {
    for (std::size_t i = 0; i < A.NumRows(); ++i)
    {
      result = 0.0;
      for (std::size_t j = 0; j < A.NumColumns(); ++j)
      {
        if (!A.IsZero(i, j))
        {
          result += A[i_block][i][j] * X[i_block][j];
        }
      }
      F(B[i_block][i], result);
    }
  }
}

template<class SparseMatrixPolicy>
void PrintMatrix(const SparseMatrixPolicy& matrix, std::size_t width)
{
  for (std::size_t i_block = 0; i_block < matrix.NumberOfBlocks(); ++i_block)
  {
    std::cout << "block: " << i_block << std::endl;
    for (std::size_t i = 0; i < matrix.NumRows(); ++i)
    {
      for (std::size_t j = 0; j < matrix.NumColumns(); ++j)
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
void testDenseMatrix()
{
  using FloatingPointType = typename MatrixPolicy::value_type;

  SparseMatrixPolicy a = SparseMatrixPolicy(
      SparseMatrixPolicy::Create(3)
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

  a[0][0][0] = 2;
  a[0][0][1] = -1;
  a[0][0][2] = -2;
  a[0][1][0] = -4;
  a[0][1][1] = 6;
  a[0][1][2] = 3;
  a[0][2][0] = -4;
  a[0][2][1] = -2;
  a[0][2][2] = 8;

  b[0][0] = 23;
  b[0][1] = 42;
  b[0][2] = 9;

  x = b;

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<MatrixPolicy>(b);
  CheckCopyToDevice<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(a, 0);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(a, 0, false);
  alu.Fill(0);
  for (std::size_t i = 0; i < a.NumRows(); ++i)
  {
    for (std::size_t j = 0; j < a.NumColumns(); ++j)
    {
      if (!a.IsZero(i, j))
      {
        alu[0][i][j] = a[0][i][j];
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
      a, b, x, [&](const FloatingPointType A, const FloatingPointType B) -> void { EXPECT_NEAR(A, B, 1.0e-5); });
}

template<class MatrixPolicy, class SparseMatrixPolicy, class LinearSolverPolicy>
void testRandomMatrix(std::size_t number_of_blocks)
{
  using FloatingPointType = typename MatrixPolicy::value_type;

  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy::Create(10).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (std::size_t i = 0; i < 10; ++i)
  {
    for (std::size_t j = 0; j < 10; ++j)
    {
      if (i == j || gen_bool())
      {
        builder = builder.WithElement(i, j);
      }
    }
  }

  SparseMatrixPolicy a(builder);
  MatrixPolicy b(number_of_blocks, 10, 0.0);
  MatrixPolicy x(number_of_blocks, 10, 0.0);

  for (std::size_t i = 0; i < 10; ++i)
  {
    for (std::size_t j = 0; j < 10; ++j)
    {
      if (!a.IsZero(i, j))
      {
        for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
        {
          a[i_block][i][j] = get_double();
        }
      }
    }
  }

  for (std::size_t i = 0; i < 10; ++i)
  {
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
    {
      b[i_block][i] = get_double();
    }
  }

  x = b;

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(a, 0);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(a, 0, false);
  alu.Fill(0);
  for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
  {
    for (std::size_t i = 0; i < a.NumRows(); ++i)
    {
      for (std::size_t j = 0; j < a.NumColumns(); ++j)
      {
        if (!a.IsZero(i, j))
        {
          alu[i_block][i][j] = a[i_block][i][j];
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
      a, b, x, [&](const FloatingPointType A, const FloatingPointType B) -> void { EXPECT_NEAR(A, B, 1.0e-6); });
}

template<class MatrixPolicy, class SparseMatrixPolicy, class LinearSolverPolicy>
void testExtremeInitialValue(std::size_t number_of_blocks, double initial_value)
{
  using FloatingPointType = typename MatrixPolicy::value_type;

  const unsigned int SEED = 12345;
  std::mt19937 generator(SEED);
  const double POINT_FIVE = 0.5;
  const double TWO = 2.0;

  auto gen_bool = std::bind(std::bernoulli_distribution(POINT_FIVE), generator);
  auto get_double = std::bind(std::lognormal_distribution<double>(-TWO, TWO), generator);
  const size_t SIZE = 30;

  auto builder = SparseMatrixPolicy::Create(SIZE).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (std::size_t i = 0; i < SIZE; ++i)
  {
    for (std::size_t j = 0; j < SIZE; ++j)
    {
      if (i == j || gen_bool())
      {
        builder = builder.WithElement(i, j);
      }
    }
  }

  SparseMatrixPolicy a(builder);
  MatrixPolicy b(number_of_blocks, SIZE, 0.0);
  MatrixPolicy x(number_of_blocks, SIZE, 0.0);

  for (std::size_t i = 0; i < SIZE; ++i)
  {
    for (std::size_t j = 0; j < SIZE; ++j)
    {
      if (!a.IsZero(i, j))
      {
        for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
        {
          a[i_block][i][j] = get_double();
        }
      }
    }
  }

  for (std::size_t i = 0; i < SIZE; ++i)
  {
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
    {
      b[i_block][i] = get_double();
    }
  }

  x = b;

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(a, initial_value);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(a, initial_value, false);
  for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
  {
    for (std::size_t i = 0; i < a.NumRows(); ++i)
    {
      for (std::size_t j = 0; j < a.NumColumns(); ++j)
      {
        if (!a.IsZero(i, j))
        {
          alu[i_block][i][j] = a[i_block][i][j];
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

  CheckResults<FloatingPointType, MatrixPolicy, SparseMatrixPolicy>(
      a, b, x, [&](const FloatingPointType A, const FloatingPointType B) -> void { EXPECT_NEAR(A, B, 2.0e-06); });
}

template<class MatrixPolicy, class SparseMatrixPolicy, class LinearSolverPolicy>
void testDiagonalMatrix(std::size_t number_of_blocks)
{
  using FloatingPointType = typename MatrixPolicy::value_type;

  auto get_double = std::bind(std::lognormal_distribution(-2.0, 4.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy::Create(6).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (std::size_t i = 0; i < 6; ++i)
  {
    builder = builder.WithElement(i, i);
  }

  SparseMatrixPolicy a(builder);
  MatrixPolicy b(number_of_blocks, 6, 0.0);
  MatrixPolicy x(number_of_blocks, 6, 0.0);

  for (std::size_t i = 0; i < 6; ++i)
  {
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
    {
      a[i_block][i][i] = get_double();
    }
  }

  for (std::size_t i = 0; i < 6; ++i)
  {
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
    {
      b[i_block][i] = get_double();
    }
  }

  x = b;

  // Only copy the data to the device when it is a CudaMatrix
  CheckCopyToDevice<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(a, 0);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(a, 0, false);
  alu.Fill(0);
  for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
  {
    for (std::size_t i = 0; i < a.NumRows(); ++i)
    {
      for (std::size_t j = 0; j < a.NumColumns(); ++j)
      {
        if (!a.IsZero(i, j))
        {
          alu[i_block][i][j] = a[i_block][i][j];
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
      a, b, x, [&](const FloatingPointType A, const FloatingPointType B) -> void { EXPECT_NEAR(A, B, 1.0e-5); });
}

template<class MatrixPolicy, class SparseMatrixPolicy>
void testMarkowitzReordering()
{
  const std::size_t ORDER = 50;
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  MatrixPolicy orig(ORDER, ORDER, 0);

  for (std::size_t i = 0; i < ORDER; ++i)
  {
    for (std::size_t j = 0; j < ORDER; ++j)
    {
      orig[i][j] = (i == j || gen_bool()) ? 1 : 0;
    }
  }

  auto reorder_map = micm::DiagonalMarkowitzReorder<MatrixPolicy>(orig);

  auto builder = SparseMatrixPolicy::Create(50);
  for (std::size_t i = 0; i < ORDER; ++i)
  {
    for (std::size_t j = 0; j < ORDER; ++j)
    {
      if (orig[i][j] != 0)
      {
        builder = builder.WithElement(i, j);
      }
    }
  }
  SparseMatrixPolicy orig_jac{ builder };

  builder = SparseMatrixPolicy::Create(50);
  for (std::size_t i = 0; i < ORDER; ++i)
  {
    for (std::size_t j = 0; j < ORDER; ++j)
    {
      if (orig[reorder_map[i]][reorder_map[j]] != 0)
      {
        builder = builder.WithElement(i, j);
      }
    }
  }
  SparseMatrixPolicy reordered_jac{ builder };

  auto orig_lu_calc = micm::LuDecompositionInPlace::Create<SparseMatrixPolicy>(orig_jac);
  auto reordered_lu_calc = micm::LuDecompositionInPlace::Create<SparseMatrixPolicy>(reordered_jac);

  auto orig_lu = orig_lu_calc.template GetLUMatrix<SparseMatrixPolicy>(orig_jac, 0.0, false);
  auto reordered_lu = reordered_lu_calc.template GetLUMatrix<SparseMatrixPolicy>(reordered_jac, 0.0, false);

  std::size_t sum_orig = 0;
  std::size_t sum_reordered = 0;
  for (std::size_t i = 0; i < reorder_map.size(); ++i)
  {
    sum_orig += i;
    sum_reordered += reorder_map[i];
  }

  EXPECT_EQ(sum_orig, sum_reordered);
  EXPECT_GT(orig_lu.FlatBlockSize(), reordered_lu.FlatBlockSize());
}
