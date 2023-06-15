#include <gtest/gtest.h>

#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/sparse_matrix.hpp>

template<class T>
void check_results(const micm::SparseMatrix<T>& A, const micm::SparseMatrix<T>& L, const micm::SparseMatrix<T>& U)
{
  EXPECT_EQ(A.size(), L.size());
  EXPECT_EQ(A.size(), U.size());
  for (std::size_t i_block = 0; i_block < A.size(); ++i_block)
  {
    for (std::size_t i = 0; i < A[i_block].size(); ++i)
    {
      for (std::size_t j = 0; j < A[i_block].size(); ++j)
      {
        T result{};
        for (std::size_t k = 0; k < A[i_block].size(); ++k)
        {
          if (!(L.IsZero(i, k) || U.IsZero(k, j)))
          {
            result += L[i_block][i][k] * U[i_block][k][j];
          }
        }
        if (A.IsZero(i, j))
        {
          EXPECT_EQ(result, 0);
        }
        else
        {
          EXPECT_EQ(A[i_block][i][j], result);
        }
      }
    }
  }
}

template<class T>
void print_matrix(const micm::SparseMatrix<T>& matrix, std::size_t width)
{
  for (std::size_t i_block = 0; i_block < matrix.size(); ++i_block)
  {
    std::cout << "block: " << i_block << std::endl;
    for (std::size_t i = 0; i < matrix[i_block].size(); ++i)
    {
      for (std::size_t j = 0; j < matrix[i_block][i].size(); ++j)
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
TEST(LuDecomposition, denseMatrix)
{
  micm::SparseMatrix<int> A = micm::SparseMatrix<int>::create(3)
                                  .with_element(0, 0)
                                  .with_element(0, 1)
                                  .with_element(0, 2)
                                  .with_element(1, 0)
                                  .with_element(1, 1)
                                  .with_element(1, 2)
                                  .with_element(2, 0)
                                  .with_element(2, 1)
                                  .with_element(2, 2);

  A[0][0][0] = 2;
  A[0][0][1] = -1;
  A[0][0][2] = -2;
  A[0][1][0] = -4;
  A[0][1][1] = 6;
  A[0][1][2] = 3;
  A[0][2][0] = -4;
  A[0][2][1] = -2;
  A[0][2][2] = 8;

  micm::LuDecomposition lud(A);
  lud.Print();

  auto LU = micm::LuDecomposition::GetLUMatrices(A);
  lud.Decompose(A, LU.first, LU.second);

  std::cout << "A matrix" << std::endl;
  print_matrix(A, 3);
  std::cout << "L matrix" << std::endl;
  print_matrix(LU.first, 3);
  std::cout << "U matrix" << std::endl;
  print_matrix(LU.second, 3);

  check_results(A, LU.first, LU.second);
}