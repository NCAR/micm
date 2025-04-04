#include <micm/solver/linear_solver.hpp>

#include <gtest/gtest.h>

#include <functional>
#include <random>

// Define the following three functions that only work for the CudaMatrix; the if constexpr statement is evalauted at
// compile-time Reference: https://www.modernescpp.com/index.php/using-requires-expression-in-c-20-as-a-standalone-feature/
template<class MatrixPolicy>
void CopyToDeviceDense(MatrixPolicy& matrix)
{
  if constexpr (requires {
                  { matrix.CopyToDevice() } -> std::same_as<void>;
                })
    matrix.CopyToDevice();
}

template<class SparseMatrixPolicy>
void CopyToDeviceSparse(SparseMatrixPolicy& matrix)
{
  if constexpr (requires {
                  { matrix.CopyToDevice() } -> std::same_as<void>;
                })
    matrix.CopyToDevice();
}

template<class MatrixPolicy>
void CopyToHostDense(MatrixPolicy& matrix)
{
  if constexpr (requires {
                  { matrix.CopyToHost() } -> std::same_as<void>;
                })
    matrix.CopyToHost();
}

template<typename T, class MatrixPolicy, class SparseMatrixPolicy>
void check_results(
    const SparseMatrixPolicy A,
    const MatrixPolicy b,
    const MatrixPolicy x,
    const std::function<void(const T, const T)> f)
{
  T result;
  EXPECT_EQ(A.NumberOfBlocks(), b.NumRows());
  EXPECT_EQ(A.NumberOfBlocks(), x.NumRows());
  for (std::size_t i_block = 0; i_block < A.NumberOfBlocks(); ++i_block)
  {
    for (std::size_t i = 0; i < A.NumRows(); ++i)
    {
      result = 0.0;
      for (std::size_t j = 0; j < A.NumColumns(); ++j)
        if (!A.IsZero(i, j))
          result += A[i_block][i][j] * x[i_block][j];
      f(b[i_block][i], result);
    }
  }
}

template<class SparseMatrixPolicy>
void print_matrix(const SparseMatrixPolicy& matrix, std::size_t width)
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
  CopyToDeviceDense<MatrixPolicy>(b);
  CopyToDeviceDense<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(A, 0);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(A, 0);
  alu.Fill(0);
  for (std::size_t i = 0; i < A.NumRows(); ++i)
    for (std::size_t j = 0; j < A.NumColumns(); ++j)
      if (!A.IsZero(i, j))
        alu[0][i][j] = A[0][i][j];

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<SparseMatrixPolicy>(alu);

  solver.Factor(alu);
  solver.template Solve<MatrixPolicy>(x, alu);

  // Only copy the data to the host when it is a CudaMatrix
  CopyToHostDense<MatrixPolicy>(x);

  check_results<FloatingPointType, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const FloatingPointType a, const FloatingPointType b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}

template<class MatrixPolicy, class SparseMatrixPolicy, class LinearSolverPolicy>
void testRandomMatrix(std::size_t number_of_blocks)
{
  using FloatingPointType = typename MatrixPolicy::value_type;

  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy::Create(10).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (i == j || gen_bool())
        builder = builder.WithElement(i, j);

  SparseMatrixPolicy A(builder);
  MatrixPolicy b(number_of_blocks, 10, 0.0);
  MatrixPolicy x(number_of_blocks, 10, 0.0);

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (!A.IsZero(i, j))
        for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
          A[i_block][i][j] = get_double();

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
      b[i_block][i] = get_double();

  x = b;

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceDense<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(A, 0);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(A, 0);
  alu.Fill(0);
  for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
    for (std::size_t i = 0; i < A.NumRows(); ++i)
      for (std::size_t j = 0; j < A.NumColumns(); ++j)
        if (!A.IsZero(i, j))
          alu[i_block][i][j] = A[i_block][i][j];

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<SparseMatrixPolicy>(alu);

  solver.Factor(alu);
  solver.template Solve<MatrixPolicy>(x, alu);

  // Only copy the data to the host when it is a CudaMatrix
  CopyToHostDense<MatrixPolicy>(x);

  check_results<FloatingPointType, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const FloatingPointType a, const FloatingPointType b) -> void { EXPECT_NEAR(a, b, 1.0e-6); });
}

template<class MatrixPolicy, class SparseMatrixPolicy, class LinearSolverPolicy>
void testExtremeInitialValue(std::size_t number_of_blocks, double initial_value)
{
  using FloatingPointType = typename MatrixPolicy::value_type;

  const unsigned int seed = 12345;
  std::mt19937 generator(seed);
  const double point_five = 0.5;
  const double two = 2.0;

  auto gen_bool = std::bind(std::bernoulli_distribution(point_five), generator);
  auto get_double = std::bind(std::lognormal_distribution<double>(-two, two), generator);
  const size_t size = 30;

  auto builder = SparseMatrixPolicy::Create(size).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (std::size_t i = 0; i < size; ++i)
  {
    for (std::size_t j = 0; j < size; ++j)
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

  for (std::size_t i = 0; i < size; ++i)
    for (std::size_t j = 0; j < size; ++j)
      if (!A.IsZero(i, j))
        for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
          A[i_block][i][j] = get_double();

  for (std::size_t i = 0; i < size; ++i)
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
      b[i_block][i] = get_double();

  x = b;

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceDense<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(A, initial_value);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(A, initial_value);
  for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
    for (std::size_t i = 0; i < A.NumRows(); ++i)
      for (std::size_t j = 0; j < A.NumColumns(); ++j)
        if (!A.IsZero(i, j))
          alu[i_block][i][j] = A[i_block][i][j];
        else
          if (!alu.IsZero(i,j)) alu[i_block][i][j] = 0;

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<SparseMatrixPolicy>(alu);

  solver.Factor(alu);

  // Only copy the data to the host when it is a CudaMatrix
  CopyToHostDense<SparseMatrixPolicy>(alu);

  solver.template Solve<MatrixPolicy>(x, alu);

  // Only copy the data to the host when it is a CudaMatrix
  CopyToHostDense<MatrixPolicy>(x);

  check_results<FloatingPointType, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const FloatingPointType a, const FloatingPointType b) -> void { EXPECT_NEAR(a, b, 2.0e-06); });
}

template<class MatrixPolicy, class SparseMatrixPolicy, class LinearSolverPolicy>
void testDiagonalMatrix(std::size_t number_of_blocks)
{
  using FloatingPointType = typename MatrixPolicy::value_type;

  auto get_double = std::bind(std::lognormal_distribution(-2.0, 4.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy::Create(6).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (std::size_t i = 0; i < 6; ++i)
    builder = builder.WithElement(i, i);

  SparseMatrixPolicy A(builder);
  MatrixPolicy b(number_of_blocks, 6, 0.0);
  MatrixPolicy x(number_of_blocks, 6, 0.0);

  for (std::size_t i = 0; i < 6; ++i)
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
      A[i_block][i][i] = get_double();

  for (std::size_t i = 0; i < 6; ++i)
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
      b[i_block][i] = get_double();

  x = b;

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceDense<MatrixPolicy>(x);

  LinearSolverPolicy solver = LinearSolverPolicy(A, 0);
  auto alu = micm::LuDecompositionInPlace::GetLUMatrix<SparseMatrixPolicy>(A, 0);
  alu.Fill(0);
  for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
    for (std::size_t i = 0; i < A.NumRows(); ++i)
      for (std::size_t j = 0; j < A.NumColumns(); ++j)
        if (!A.IsZero(i, j))
          alu[i_block][i][j] = A[i_block][i][j];

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<SparseMatrixPolicy>(alu);

  solver.Factor(alu);
  solver.template Solve<MatrixPolicy>(x, alu);

  // Only copy the data to the host when it is a CudaMatrix
  CopyToHostDense<MatrixPolicy>(x);

  check_results<FloatingPointType, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const FloatingPointType a, const FloatingPointType b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}

template<class MatrixPolicy, class SparseMatrixPolicy>
void testMarkowitzReordering()
{
  const std::size_t order = 50;
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  MatrixPolicy orig(order, order, 0);

  for (std::size_t i = 0; i < order; ++i)
    for (std::size_t j = 0; j < order; ++j)
      orig[i][j] = (i == j || gen_bool()) ? 1 : 0;

  auto reorder_map = micm::DiagonalMarkowitzReorder<MatrixPolicy>(orig);

  auto builder = SparseMatrixPolicy::Create(50);
  for (std::size_t i = 0; i < order; ++i)
    for (std::size_t j = 0; j < order; ++j)
      if (orig[i][j] != 0)
        builder = builder.WithElement(i, j);
  SparseMatrixPolicy orig_jac{ builder };

  builder = SparseMatrixPolicy::Create(50);
  for (std::size_t i = 0; i < order; ++i)
    for (std::size_t j = 0; j < order; ++j)
      if (orig[reorder_map[i]][reorder_map[j]] != 0)
        builder = builder.WithElement(i, j);
  SparseMatrixPolicy reordered_jac{ builder };

  auto orig_LU_calc = micm::LuDecompositionInPlace::Create<SparseMatrixPolicy>(orig_jac);
  auto reordered_LU_calc = micm::LuDecompositionInPlace::Create<SparseMatrixPolicy>(reordered_jac);

  auto orig_LU = orig_LU_calc.template GetLUMatrix<SparseMatrixPolicy>(orig_jac, 0.0);
  auto reordered_LU = reordered_LU_calc.template GetLUMatrix<SparseMatrixPolicy>(reordered_jac, 0.0);

  std::size_t sum_orig = 0;
  std::size_t sum_reordered = 0;
  for (std::size_t i = 0; i < reorder_map.size(); ++i)
  {
    sum_orig += i;
    sum_reordered += reorder_map[i];
  }

  EXPECT_EQ(sum_orig, sum_reordered);
  EXPECT_GT(orig_LU.FlatBlockSize(), reordered_LU.FlatBlockSize());
}
