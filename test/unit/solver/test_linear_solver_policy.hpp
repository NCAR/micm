#include <micm/solver/linear_solver.hpp>
#include <micm/solver/lu_decomposition.hpp>

#include <gtest/gtest.h>

#include <functional>
#include <random>

// Define the following three functions that only work for the CudaMatrix; the if constexpr statement is evalauted at
// compile-time Reference: https://www.modernescpp.com/index.php/using-requires-expression-in-c-20-as-a-standalone-feature/
template<typename T, template<class> class MatrixPolicy>
void CopyToDeviceDense(MatrixPolicy<T>& matrix)
{
  if constexpr (requires {
                  {
                    matrix.CopyToDevice()
                    } -> std::same_as<void>;
                })
    matrix.CopyToDevice();
}

template<typename T, template<class> class SparseMatrixPolicy>
void CopyToDeviceSparse(SparseMatrixPolicy<T>& matrix)
{
  if constexpr (requires {
                  {
                    matrix.CopyToDevice()
                    } -> std::same_as<void>;
                })
    matrix.CopyToDevice();
}

template<typename T, template<class> class MatrixPolicy>
void CopyToHostDense(MatrixPolicy<T>& matrix)
{
  if constexpr (requires {
                  {
                    matrix.CopyToHost()
                    } -> std::same_as<void>;
                })
    matrix.CopyToHost();
}

template<typename T, template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
void check_results(
    const SparseMatrixPolicy<T> A,
    const MatrixPolicy<T> b,
    const MatrixPolicy<T> x,
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

template<typename T, template<class> class SparseMatrixPolicy>
void print_matrix(const SparseMatrixPolicy<T>& matrix, std::size_t width)
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
          std::cout << " " << std::setfill('-') << std::setw(width) << "-"
                    << " ";
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

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
void testDenseMatrix(const std::function<LinearSolverPolicy(const SparseMatrixPolicy<double>, double)> create_linear_solver)
{
  SparseMatrixPolicy<double> A = SparseMatrixPolicy<double>(SparseMatrixPolicy<double>::Create(3)
                                                                .InitialValue(1.0e-30)
                                                                .WithElement(0, 0)
                                                                .WithElement(0, 1)
                                                                .WithElement(0, 2)
                                                                .WithElement(1, 0)
                                                                .WithElement(1, 1)
                                                                .WithElement(1, 2)
                                                                .WithElement(2, 0)
                                                                .WithElement(2, 1)
                                                                .WithElement(2, 2));
  MatrixPolicy<double> b(1, 3, 0.0);
  MatrixPolicy<double> x(1, 3, 100.0);

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

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<double, SparseMatrixPolicy>(A);
  CopyToDeviceDense<double, MatrixPolicy>(b);
  CopyToDeviceDense<double, MatrixPolicy>(x);

  LinearSolverPolicy solver = create_linear_solver(A, 1.0e-30);
  auto lu = micm::LuDecomposition::GetLUMatrices<double, SparseMatrixPolicy>(A, 1.0e-30);
  auto lower_matrix = std::move(lu.first);
  auto upper_matrix = std::move(lu.second);

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<double, SparseMatrixPolicy>(lower_matrix);
  CopyToDeviceSparse<double, SparseMatrixPolicy>(upper_matrix);

  solver.Factor(A, lower_matrix, upper_matrix);
  solver.template Solve<MatrixPolicy>(b, x, lower_matrix, upper_matrix);

  // Only copy the data to the host when it is a CudaMatrix
  CopyToHostDense<double, MatrixPolicy>(x);

  check_results<double, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
void testRandomMatrix(
    const std::function<LinearSolverPolicy(const SparseMatrixPolicy<double>, double)> create_linear_solver,
    std::size_t number_of_blocks)
{
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy<double>::Create(10).SetNumberOfBlocks(number_of_blocks).InitialValue(1.0e-30);
  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (i == j || gen_bool())
        builder = builder.WithElement(i, j);

  SparseMatrixPolicy<double> A(builder);
  MatrixPolicy<double> b(number_of_blocks, 10, 0.0);
  MatrixPolicy<double> x(number_of_blocks, 10, 100.0);

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (!A.IsZero(i, j))
        for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
          A[i_block][i][j] = get_double();

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
      b[i_block][i] = get_double();

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<double, SparseMatrixPolicy>(A);
  CopyToDeviceDense<double, MatrixPolicy>(b);
  CopyToDeviceDense<double, MatrixPolicy>(x);

  LinearSolverPolicy solver = create_linear_solver(A, 1.0e-30);
  auto lu = micm::LuDecomposition::GetLUMatrices<double, SparseMatrixPolicy>(A, 1.0e-30);
  auto lower_matrix = std::move(lu.first);
  auto upper_matrix = std::move(lu.second);

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<double, SparseMatrixPolicy>(lower_matrix);
  CopyToDeviceSparse<double, SparseMatrixPolicy>(upper_matrix);

  solver.Factor(A, lower_matrix, upper_matrix);
  solver.template Solve<MatrixPolicy>(b, x, lower_matrix, upper_matrix);

  // Only copy the data to the host when it is a CudaMatrix
  CopyToHostDense<double, MatrixPolicy>(x);

  check_results<double, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
void testDiagonalMatrix(
    const std::function<LinearSolverPolicy(const SparseMatrixPolicy<double>, double)> create_linear_solver,
    std::size_t number_of_blocks)
{
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 4.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy<double>::Create(6).SetNumberOfBlocks(number_of_blocks).InitialValue(1.0e-30);
  for (std::size_t i = 0; i < 6; ++i)
    builder = builder.WithElement(i, i);

  SparseMatrixPolicy<double> A(builder);
  MatrixPolicy<double> b(number_of_blocks, 6, 0.0);
  MatrixPolicy<double> x(number_of_blocks, 6, 100.0);

  for (std::size_t i = 0; i < 6; ++i)
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
      A[i_block][i][i] = get_double();

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<double, SparseMatrixPolicy>(A);
  CopyToDeviceDense<double, MatrixPolicy>(b);
  CopyToDeviceDense<double, MatrixPolicy>(x);

  LinearSolverPolicy solver = create_linear_solver(A, 1.0e-30);
  auto lu = micm::LuDecomposition::GetLUMatrices<double, SparseMatrixPolicy>(A, 1.0e-30);
  auto lower_matrix = std::move(lu.first);
  auto upper_matrix = std::move(lu.second);

  // Only copy the data to the device when it is a CudaMatrix
  CopyToDeviceSparse<double, SparseMatrixPolicy>(lower_matrix);
  CopyToDeviceSparse<double, SparseMatrixPolicy>(upper_matrix);

  solver.Factor(A, lower_matrix, upper_matrix);
  solver.template Solve<MatrixPolicy>(b, x, lower_matrix, upper_matrix);

  // Only copy the data to the host when it is a CudaMatrix
  CopyToHostDense<double, MatrixPolicy>(x);

  check_results<double, MatrixPolicy, SparseMatrixPolicy>(
      A, b, x, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-5); });
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
void testMarkowitzReordering()
{
  const std::size_t order = 50;
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  MatrixPolicy<int> orig(order, order, 0);

  for (std::size_t i = 0; i < order; ++i)
    for (std::size_t j = 0; j < order; ++j)
      orig[i][j] = (i == j || gen_bool()) ? 1 : 0;

  auto reorder_map = micm::DiagonalMarkowitzReorder<MatrixPolicy>(orig);

  auto builder = SparseMatrixPolicy<double>::Create(50);
  for (std::size_t i = 0; i < order; ++i)
    for (std::size_t j = 0; j < order; ++j)
      if (orig[i][j] != 0)
        builder = builder.WithElement(i, j);
  SparseMatrixPolicy<double> orig_jac{ builder };

  builder = SparseMatrixPolicy<double>::Create(50);
  for (std::size_t i = 0; i < order; ++i)
    for (std::size_t j = 0; j < order; ++j)
      if (orig[reorder_map[i]][reorder_map[j]] != 0)
        builder = builder.WithElement(i, j);
  SparseMatrixPolicy<double> reordered_jac{ builder };

  auto orig_LU_calc = micm::LuDecomposition::Create<double, SparseMatrixPolicy>(orig_jac);
  auto reordered_LU_calc = micm::LuDecomposition::Create<double, SparseMatrixPolicy>(reordered_jac);

  auto orig_LU = orig_LU_calc.template GetLUMatrices<double, SparseMatrixPolicy>(orig_jac, 0.0);
  auto reordered_LU = reordered_LU_calc.template GetLUMatrices<double, SparseMatrixPolicy>(reordered_jac, 0.0);

  std::size_t sum_orig = 0;
  std::size_t sum_reordered = 0;
  for (std::size_t i = 0; i < reorder_map.size(); ++i)
  {
    sum_orig += i;
    sum_reordered += reorder_map[i];
  }

  EXPECT_EQ(sum_orig, sum_reordered);
  EXPECT_GT(
      orig_LU.first.RowIdsVector().size() + orig_LU.second.RowIdsVector().size(),
      reordered_LU.first.RowIdsVector().size() + reordered_LU.second.RowIdsVector().size());
}