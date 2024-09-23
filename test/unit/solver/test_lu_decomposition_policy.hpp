#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <gtest/gtest.h>

#include <functional>
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
  std::cout << "A.NumRows: " << A.NumRows() << " A.NumColumns: " A.NumColumns() << std::endl;
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
          // f(result, T{});
        }
        else
        {
          f(result, A[i_block][i][j]);
        }
      }
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

// tests example from https://www.geeksforgeeks.org/doolittle-algorithm-lu-decomposition/
template<class SparseMatrixPolicy, class LuDecompositionPolicy>
void testDenseMatrix()
{
  SparseMatrixPolicy A = SparseMatrixPolicy(SparseMatrixPolicy::Create(3)
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

  A[0][0][0] = 2;
  A[0][0][1] = -1;
  A[0][0][2] = -2;
  A[0][1][0] = -4;
  A[0][1][1] = 6;
  A[0][1][2] = 3;
  A[0][2][0] = -4;
  A[0][2][1] = -2;
  A[0][2][2] = 8;

  LuDecompositionPolicy lud = LuDecompositionPolicy(A);
  auto LU = micm::LuDecomposition::GetLUMatrices<SparseMatrixPolicy>(A, 1.0e-30);
  bool is_singular{ false };
  lud.template Decompose<SparseMatrixPolicy>(A, LU.first, LU.second, is_singular);
  check_results<double, SparseMatrixPolicy>(
      A, LU.first, LU.second, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-10); });
}

template<class SparseMatrixPolicy, class LuDecompositionPolicy>
void testSingularMatrix()
{
  SparseMatrixPolicy A = SparseMatrixPolicy(
      SparseMatrixPolicy::Create(2).InitialValue(1.0e-30).WithElement(0, 0).WithElement(0, 1).WithElement(1, 0).WithElement(
          1, 1));

  A[0][0][0] = 0;
  A[0][0][1] = 1;
  A[0][1][0] = 1;
  A[0][1][1] = 1;

  LuDecompositionPolicy lud = LuDecompositionPolicy(A);
  auto LU = micm::LuDecomposition::GetLUMatrices<SparseMatrixPolicy>(A, 1.0E-30);
  bool is_singular{ false };
  lud.template Decompose<SparseMatrixPolicy>(A, LU.first, LU.second, is_singular);
  EXPECT_TRUE(is_singular);
  A[0][0][0] = 12;
  lud.template Decompose<SparseMatrixPolicy>(A, LU.first, LU.second, is_singular);
  EXPECT_FALSE(is_singular);
}

template<class SparseMatrixPolicy, class LuDecompositionPolicy>
void testRandomMatrix(std::size_t number_of_blocks)
{
  auto gen_bool = std::bind(std::uniform_int_distribution<>(0, 1), std::default_random_engine());
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy::Create(10).SetNumberOfBlocks(number_of_blocks).InitialValue(1.0e-30);
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

  LuDecompositionPolicy lud = LuDecompositionPolicy(A);
  auto LU = micm::LuDecomposition::GetLUMatrices<SparseMatrixPolicy>(A, 1.0e-30);
  bool is_singular{ false };
  lud.template Decompose<SparseMatrixPolicy>(A, LU.first, LU.second, is_singular);
  check_results<double, SparseMatrixPolicy>(
      A, LU.first, LU.second, [&](const double a, const double b) -> void { 
        EXPECT_LT(std::abs((a-b)/b), 1.0e-10); 
      });
}

template<class SparseMatrixPolicy, class LuDecompositionPolicy>
void testDiagonalMatrix(std::size_t number_of_blocks)
{
  auto get_double = std::bind(std::lognormal_distribution(-2.0, 4.0), std::default_random_engine());

  auto builder = SparseMatrixPolicy::Create(6).SetNumberOfBlocks(number_of_blocks).InitialValue(1.0e-30);
  for (std::size_t i = 0; i < 6; ++i)
    builder = builder.WithElement(i, i);

  SparseMatrixPolicy A(builder);

  for (std::size_t i = 0; i < 6; ++i)
    for (std::size_t i_block = 0; i_block < number_of_blocks; ++i_block)
      A[i_block][i][i] = get_double();

  LuDecompositionPolicy lud = LuDecompositionPolicy(A);
  auto LU = micm::LuDecomposition::GetLUMatrices<SparseMatrixPolicy>(A, 1.0e-30);
  bool is_singular{ false };
  lud.template Decompose<SparseMatrixPolicy>(A, LU.first, LU.second, is_singular);
  check_results<double, SparseMatrixPolicy>(
      A, LU.first, LU.second, [&](const double a, const double b) -> void { EXPECT_NEAR(a, b, 1.0e-10); });
}