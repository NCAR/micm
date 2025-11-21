#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

#include <functional>
#include <iomanip>
#include <random>

template<typename T, class SparseMatrixPolicy>
void check_results(
    const SparseMatrixPolicy& A,
    const SparseMatrixPolicy& L,
    const SparseMatrixPolicy& U,
    const std::function<void(const T, const T)> f)
{
  EXPECT_EQ(A.NumberOfBlocks(), L.NumberOfBlocks());
  EXPECT_EQ(A.NumberOfBlocks(), U.NumberOfBlocks());
  for (std::size_t i_block = 0; i_block < A.NumberOfBlocks(); ++i_block)
  {
    for (std::size_t i = 0; i < A.NumRows(); ++i)
    {
      for (std::size_t j = 0; j < A.NumColumns(); ++j)
      {
        T result{};
        for (std::size_t k = 0; k < A.NumRows(); ++k)
        {
          if (!(L.IsZero(i, k) || U.IsZero(k, j)))
          {
            result += L[i_block][i][k] * U[i_block][k][j];
          }
        }
        // Make sure these are actually triangular matrices
        EXPECT_TRUE(i >= j || L.IsZero(i, j));
        EXPECT_TRUE(j >= i || U.IsZero(i, j));
        if (A.IsZero(i, j))
        {
          f(result, T{});
        }
        else
        {
          f(result, A[i_block][i][j]);
        }
      }
    }
  }
}

void print_matrix(const auto& matrix, std::size_t width)
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

// tests example from https://www.geeksforgeeks.org/doolittle-algorithm-lu-decomposition/
template<class SparseMatrixPolicy, class LuDecompositionPolicy>
void testDenseMatrix()
{
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

  A[0][0][0] = 2;
  A[0][0][1] = -1;
  A[0][0][2] = -2;
  A[0][1][0] = -4;
  A[0][1][1] = 6;
  A[0][1][2] = 3;
  A[0][2][0] = -4;
  A[0][2][1] = -2;
  A[0][2][2] = 8;

  LuDecompositionPolicy lud =
      LuDecompositionPolicy::template Create<SparseMatrixPolicy, SparseMatrixPolicy, SparseMatrixPolicy>(A);
  auto LU =
      LuDecompositionPolicy::template GetLUMatrices<SparseMatrixPolicy, SparseMatrixPolicy, SparseMatrixPolicy>(A, 0, false);
  lud.template Decompose<SparseMatrixPolicy>(A, LU.first, LU.second);
  check_results<double, SparseMatrixPolicy>(
      A, LU.first, LU.second, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-10); });
}

template<class SparseMatrixPolicy, class LuDecompositionPolicy>
void testRandomMatrix(std::size_t number_of_blocks)
{
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy::Create(10).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (i == j || gen_bool())
        builder = builder.WithElement(i, j);

  SparseMatrixPolicy A(builder);

  for (std::size_t i = 0; i < 10; ++i)
    for (std::size_t j = 0; j < 10; ++j)
      if (!A.IsZero(i, j))
        for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
          A[i_block][i][j] = get_double();

  LuDecompositionPolicy lud =
      LuDecompositionPolicy::template Create<SparseMatrixPolicy, SparseMatrixPolicy, SparseMatrixPolicy>(A);
  auto LU =
      LuDecompositionPolicy::template GetLUMatrices<SparseMatrixPolicy, SparseMatrixPolicy, SparseMatrixPolicy>(A, 0, false);
  lud.template Decompose<SparseMatrixPolicy>(A, LU.first, LU.second);
  check_results<double, SparseMatrixPolicy>(
      A, LU.first, LU.second, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-9); });
}

template<class SparseMatrixPolicy, class LuDecompositionPolicy>
void testExtremeValueInitialization(std::size_t number_of_blocks, double initial_value)
{
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());
  auto size = 10;

  auto builder = SparseMatrixPolicy::Create(10).SetNumberOfBlocks(number_of_blocks).InitialValue(initial_value);
  for (std::size_t i = 0; i < size; ++i)
    for (std::size_t j = 0; j < size; ++j)
      if (i == j || gen_bool())
        builder = builder.WithElement(i, j);

  SparseMatrixPolicy A(builder);

  // for nvhpc, the lognormal distribution produces significantly different values
  // for very large numbers of grid cells
  // To keep the accuracy on the check results function small, we only generat 1 blocks worth of
  // random values and then copy that into every other block
  for (std::size_t i = 0; i < size; ++i)
    for (std::size_t j = 0; j < size; ++j)
      if (!A.IsZero(i, j))
      {
        A[0][i][j] = get_double();
        for (std::size_t i_block = 1; i_block < number_of_blocks; ++i_block)
          A[i_block][i][j] = A[0][i][j];
      }

  LuDecompositionPolicy lud =
      LuDecompositionPolicy::template Create<SparseMatrixPolicy, SparseMatrixPolicy, SparseMatrixPolicy>(A);

  auto LU = LuDecompositionPolicy::template GetLUMatrices<SparseMatrixPolicy, SparseMatrixPolicy, SparseMatrixPolicy>(
      A, initial_value, false);

  CheckCopyToDevice<SparseMatrixPolicy>(A);
  CheckCopyToDevice<SparseMatrixPolicy>(LU.first);
  CheckCopyToDevice<SparseMatrixPolicy>(LU.second);

  lud.template Decompose<SparseMatrixPolicy>(A, LU.first, LU.second);

  CheckCopyToHost<SparseMatrixPolicy>(LU.first);
  CheckCopyToHost<SparseMatrixPolicy>(LU.second);

  check_results<double, SparseMatrixPolicy>(
      A, LU.first, LU.second, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-09); });
}

template<class SparseMatrixPolicy, class LuDecompositionPolicy>
void testDiagonalMatrix(std::size_t number_of_blocks)
{
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 4.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy::Create(6).SetNumberOfBlocks(number_of_blocks).InitialValue(0);
  for (std::size_t i = 0; i < 6; ++i)
    builder = builder.WithElement(i, i);

  SparseMatrixPolicy A(builder);

  for (std::size_t i = 0; i < 6; ++i)
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
      A[i_block][i][i] = get_double();

  LuDecompositionPolicy lud =
      LuDecompositionPolicy::template Create<SparseMatrixPolicy, SparseMatrixPolicy, SparseMatrixPolicy>(A);
  auto LU =
      LuDecompositionPolicy::template GetLUMatrices<SparseMatrixPolicy, SparseMatrixPolicy, SparseMatrixPolicy>(A, 0, false);
  lud.template Decompose<SparseMatrixPolicy>(A, LU.first, LU.second);
  check_results<double, SparseMatrixPolicy>(
      A, LU.first, LU.second, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-10); });
}
