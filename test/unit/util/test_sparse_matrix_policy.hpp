#include <micm/util/sparse_matrix.hpp>

#include <gtest/gtest.h>

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<double, OrderingPolicy> testZeroMatrix()
{
  auto builder = MatrixPolicy<double, OrderingPolicy>::Create(3);

  EXPECT_EQ(builder.NumberOfElements(), 0);

  MatrixPolicy<double, OrderingPolicy> matrix{ builder };

  EXPECT_EQ(matrix.FlatBlockSize(), 0);

  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(0, 0); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ZeroElementAccess));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(6, 0); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 3); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(6, 3); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { bool isZero = matrix.IsZero(6, 0); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { bool isZero = matrix.IsZero(1, 3); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { bool isZero = matrix.IsZero(6, 3); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[0][0][4] = 2.0; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[1][0][0] = 2.0; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[0][3][0] = 2.0; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[0][1][1] = 2.0; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ZeroElementAccess));
        throw;
      },
      std::system_error);
  return matrix;
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<double, OrderingPolicy> testConstZeroMatrix()
{
  auto builder = MatrixPolicy<double, OrderingPolicy>::Create(3);

  EXPECT_EQ(builder.NumberOfElements(), 0);

  const MatrixPolicy<double, OrderingPolicy> matrix{ builder };

  EXPECT_EQ(matrix.FlatBlockSize(), 0);

  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(0, 0); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ZeroElementAccess));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(6, 0); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 3); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(6, 3); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { bool isZero = matrix.IsZero(6, 0); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { bool isZero = matrix.IsZero(1, 3); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { bool isZero = matrix.IsZero(6, 3); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  return matrix;
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<double, OrderingPolicy> testSingleBlockMatrix()
{
  auto builder = MatrixPolicy<double, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(3, 2)
                     .WithElement(2, 3)
                     .WithElement(1, 1)
                     .WithElement(2, 1);
  // 0 X 0 0
  // 0 X 0 0
  // 0 X 0 X
  // 0 0 X 0
  {
    EXPECT_EQ(builder.NumberOfElements(), 5);
  }

  MatrixPolicy<double, OrderingPolicy> matrix{ builder };

  {
    auto diagonal_ids = matrix.DiagonalIndices(0);
    EXPECT_EQ(diagonal_ids.size(), 1);
    EXPECT_EQ(diagonal_ids[0], matrix.VectorIndex(1, 1));
  }

  EXPECT_EQ(matrix.FlatBlockSize(), 5);

  EXPECT_EQ(matrix.IsZero(0, 1), false);
  EXPECT_EQ(matrix.IsZero(3, 2), false);
  EXPECT_EQ(matrix.IsZero(2, 1), false);
  EXPECT_EQ(matrix.IsZero(2, 2), true);
  EXPECT_EQ(matrix.IsZero(0, 0), true);
  EXPECT_EQ(matrix.IsZero(3, 3), true);

  EXPECT_EQ(matrix[0][2][1], 0);
  matrix[0][2][1] = 45;
  EXPECT_EQ(matrix[0][2][1], 45);

  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(4, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 5); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 0, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(2, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ZeroElementAccess));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(0, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ZeroElementAccess));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[0][0][4] = 2; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[1][0][0] = 2; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[0][5][0] = 2; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[0][3][3] = 2; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ZeroElementAccess));
        throw;
      },
      std::system_error);
  return matrix;
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<int, OrderingPolicy> testConstSingleBlockMatrix()
{
  auto builder = MatrixPolicy<int, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(3, 2)
                     .WithElement(2, 3)
                     .WithElement(1, 1)
                     .WithElement(2, 1);
  // 0 X 0 0
  // 0 X 0 0
  // 0 X 0 X
  // 0 0 X 0
  MatrixPolicy<int, OrderingPolicy> orig_matrix{ builder };
  orig_matrix[0][2][1] = 45;
  orig_matrix[0][3][2] = 42;
  orig_matrix[0][2][3] = 21;
  const MatrixPolicy<int, OrderingPolicy> matrix = orig_matrix;

  {
    auto diagonal_ids = matrix.DiagonalIndices(0);
    EXPECT_EQ(diagonal_ids.size(), 1);
    EXPECT_EQ(diagonal_ids[0], matrix.VectorIndex(1, 1));
  }

  EXPECT_EQ(matrix.FlatBlockSize(), 5);

  EXPECT_EQ(matrix.IsZero(0, 1), false);
  EXPECT_EQ(matrix.IsZero(3, 2), false);
  EXPECT_EQ(matrix.IsZero(2, 1), false);
  EXPECT_EQ(matrix.IsZero(2, 2), true);
  EXPECT_EQ(matrix.IsZero(0, 0), true);
  EXPECT_EQ(matrix.IsZero(3, 3), true);

  EXPECT_EQ(matrix[0][2][1], 45);

  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(4, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 5); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 0, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(2, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ZeroElementAccess));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(0, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ZeroElementAccess));
        throw;
      },
      std::system_error);
  return matrix;
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<double, OrderingPolicy> testMultiBlockMatrix()
{
  auto builder = MatrixPolicy<double, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(3, 2)
                     .WithElement(2, 3)
                     .WithElement(1, 1)
                     .WithElement(2, 1)
                     .InitialValue(24)
                     .SetNumberOfBlocks(53);
  // 0 X 0 0
  // 0 X 0 0
  // 0 X 0 X
  // 0 0 X 0

  EXPECT_EQ(builder.NumberOfElements(), 5 * 53);

  MatrixPolicy<double, OrderingPolicy> matrix{ builder };

  EXPECT_EQ(matrix.FlatBlockSize(), 5);

  EXPECT_EQ(matrix[0][2][1], 24);
  matrix[0][2][1] = 45;
  EXPECT_EQ(matrix[0][2][1], 45);

  matrix[43][3][2] = 29;
  EXPECT_EQ(matrix[43][3][2], 29);

  matrix[33][2][3] = 79;
  EXPECT_EQ(matrix[33][2][3], 79);
  
  matrix[2][2][3] = 64;
  EXPECT_EQ(matrix[2][2][3], 64);

  auto diagonal_ids = matrix.DiagonalIndices(0);
  EXPECT_EQ(diagonal_ids.size(), 1);
  EXPECT_EQ(diagonal_ids[0], matrix.VectorIndex(0, 1, 1));
  diagonal_ids = matrix.DiagonalIndices(1);
  EXPECT_EQ(diagonal_ids.size(), 1);
  EXPECT_EQ(diagonal_ids[0], matrix.VectorIndex(1, 1, 1));
  diagonal_ids = matrix.DiagonalIndices(2);
  EXPECT_EQ(diagonal_ids.size(), 1);
  EXPECT_EQ(diagonal_ids[0], matrix.VectorIndex(2, 1, 1));

  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(0, 4, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(2, 1, 5); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(54, 0, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 2, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ZeroElementAccess));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(2, 0, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ZeroElementAccess));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(0, 1); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::MissingBlockIndex));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[0][0][4] = 2; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[53][0][0] = 2; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[0][5][0] = 2; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { matrix[0][3][3] = 2; } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ZeroElementAccess));
        throw;
      },
      std::system_error);
  return matrix;
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<double, OrderingPolicy> testSetScalar()
{
  auto builder = MatrixPolicy<double, OrderingPolicy>::Create(3);

  MatrixPolicy<double, OrderingPolicy> matrix{ builder };

  matrix = 2.0;

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 2.0);
  }

  return matrix;
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<int, OrderingPolicy> testAddToDiagonal()
{
  auto builder = MatrixPolicy<int, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(3, 2)
                     .WithElement(0, 1)
                     .WithElement(2, 3)
                     .WithElement(1, 1)
                     .WithElement(2, 1)
                     .InitialValue(24)
                     .SetNumberOfBlocks(3);
  // 0 X 0 0
  // 0 X 0 0
  // 0 X 0 X
  // 0 0 X 0
  MatrixPolicy<int, OrderingPolicy> matrix{ builder };

  matrix[0][1][1] = 1;
  matrix[1][1][1] = 3;
  matrix[2][1][1] = 5;

  matrix.AddToDiagonal(7);

  EXPECT_EQ(matrix[0][1][1], 8);
  EXPECT_EQ(matrix[1][1][1], 10);
  EXPECT_EQ(matrix[2][1][1], 12);
  EXPECT_EQ(matrix[0][0][1], 24);
  EXPECT_EQ(matrix[0][2][1], 24);
  EXPECT_EQ(matrix[0][2][3], 24);
  EXPECT_EQ(matrix[0][3][2], 24);
  EXPECT_EQ(matrix[1][0][1], 24);
  EXPECT_EQ(matrix[1][2][1], 24);
  EXPECT_EQ(matrix[1][2][3], 24);
  EXPECT_EQ(matrix[1][3][2], 24);
  EXPECT_EQ(matrix[2][0][1], 24);
  EXPECT_EQ(matrix[2][2][1], 24);
  EXPECT_EQ(matrix[2][2][3], 24);
  EXPECT_EQ(matrix[2][3][2], 24);

  return matrix;
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
void testPrintNonZero()
{
  auto builder = MatrixPolicy<int, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(1, 1)
                     .WithElement(2, 1)
                     .WithElement(2, 3)
                     .WithElement(3, 2)
                     .InitialValue(32)
                     .SetNumberOfBlocks(3);
  // 0 X 0 0
  // 0 X 0 0
  // 0 X 0 X
  // 0 0 X 0
  MatrixPolicy<int, OrderingPolicy> matrix{ builder };

  matrix[0][1][1] = 2;
  matrix[1][1][1] = 4;
  matrix[2][1][1] = 6;

  matrix.AddToDiagonal(5);

  std::stringstream ss;
  matrix.PrintNonZeroElements(ss);

  std::string expected_output =
      "Block 0\n"
      "0, 1, 32\n"
      "1, 1, 7\n"
      "2, 1, 32\n"
      "2, 3, 32\n"
      "3, 2, 32\n"
      "Block 1\n"
      "0, 1, 32\n"
      "1, 1, 9\n"
      "2, 1, 32\n"
      "2, 3, 32\n"
      "3, 2, 32\n"
      "Block 2\n"
      "0, 1, 32\n"
      "1, 1, 11\n"
      "2, 1, 32\n"
      "2, 3, 32\n"
      "3, 2, 32\n";

  EXPECT_EQ(ss.str(), expected_output);
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<int, OrderingPolicy> testPrint()
{
  auto builder = MatrixPolicy<int, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(3, 2)
                     .WithElement(0, 1)
                     .WithElement(2, 3)
                     .WithElement(1, 1)
                     .WithElement(2, 1)
                     .InitialValue(24)
                     .SetNumberOfBlocks(3);
  // 0 X 0 0
  // 0 X 0 0
  // 0 X 0 X
  // 0 0 X 0
  MatrixPolicy<int, OrderingPolicy> matrix{ builder };

  matrix[0][1][1] = 1;
  matrix[1][1][1] = 3;
  matrix[2][1][1] = 5;

  matrix.AddToDiagonal(7);

  std::stringstream ss, endline;
  ss << matrix;
  endline << std::endl;

  std::string expected_output = "Block 0" + endline.str() + "0,24,0,0" + endline.str() + "0,8,0,0" + endline.str() +
                                "0,24,0,24" + endline.str() + "0,0,24,0" + endline.str() + "Block 1" + endline.str() +
                                "0,24,0,0" + endline.str() + "0,10,0,0" + endline.str() + "0,24,0,24" + endline.str() +
                                "0,0,24,0" + endline.str() + "Block 2" + endline.str() + "0,24,0,0" + endline.str() +
                                "0,12,0,0" + endline.str() + "0,24,0,24" + endline.str() + "0,0,24,0" + endline.str();
  EXPECT_EQ(ss.str(), expected_output);

  return matrix;
}