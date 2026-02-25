#include <micm/util/sparse_matrix.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/vector_matrix.hpp>

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

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<double, OrderingPolicy> testArrayFunction()
{
  auto builder = MatrixPolicy<double, OrderingPolicy>::Create(4)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .WithElement(2, 2)
                     .WithElement(2, 3)
                     .WithElement(3, 3)
                     .SetNumberOfBlocks(3);
  // X 0 0 0
  // 0 X 0 0
  // 0 0 X X
  // 0 0 0 X
  MatrixPolicy<double, OrderingPolicy> matrix{ builder };

  // set some values, with unique values in different blocks
  matrix = 1.0;
  matrix[0][0][0] = 1.0;
  matrix[1][1][1] = 2.0;
  matrix[2][2][2] = 3.0;
  matrix[2][2][3] = 4.0;
  matrix[2][3][3] = 5.0;

  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
      [](auto&& mat)
      {
        auto tmp = mat.GetBlockVariable();
        mat.ForEachBlock(
            [&tmp](const double& a, const double& b, const double& c, const double& d, double& t)
            { t = a + b + c + d; },
            mat.GetConstBlockView(0,0),
            mat.GetConstBlockView(1,1),
            mat.GetConstBlockView(2,2),
            mat.GetConstBlockView(2,3),
            tmp);
        mat.ForEachBlock(
            [&tmp](double& d, const double& t)
            { d = 2.0 * t; },
            mat.GetBlockView(2,3),
            tmp);
      }, matrix); // pass matrix so the type and dimensions are known by the function

  func(matrix);

  // Check results
  EXPECT_EQ(matrix[0][2][3], 2.0 * (1.0 + 1.0 + 1.0 + 1.0)); // 8.0
  EXPECT_EQ(matrix[1][2][3], 2.0 * (1.0 + 2.0 + 1.0 + 1.0)); // 10.0
  EXPECT_EQ(matrix[2][2][3], 2.0 * (1.0 + 1.0 + 3.0 + 4.0)); // 18.0 (was 16.0 in comment, fixed)
  EXPECT_EQ(matrix[0][0][0], 1.0);
  EXPECT_EQ(matrix[1][0][0], 1.0);
  EXPECT_EQ(matrix[2][0][0], 1.0);
  EXPECT_EQ(matrix[0][1][1], 1.0);
  EXPECT_EQ(matrix[1][1][1], 2.0);
  EXPECT_EQ(matrix[2][1][1], 1.0);
  EXPECT_EQ(matrix[0][2][2], 1.0);
  EXPECT_EQ(matrix[1][2][2], 1.0);
  EXPECT_EQ(matrix[2][2][2], 3.0);
  EXPECT_EQ(matrix[0][3][3], 1.0);
  EXPECT_EQ(matrix[1][3][3], 1.0);
  EXPECT_EQ(matrix[2][3][3], 5.0);

  // Use a different matrix with the same dimensions
  MatrixPolicy<double, OrderingPolicy> matrix2{ builder };
  matrix2 = -1.0;
  func(matrix2);
  EXPECT_EQ(matrix2[0][2][3], 2.0 * (-1.0 + -1.0 + -1.0 + -1.0)); // -8.0
  EXPECT_EQ(matrix2[1][2][3], 2.0 * (-1.0 + -1.0 + -1.0 + -1.0)); // -8.0
  EXPECT_EQ(matrix2[2][2][3], 2.0 * (-1.0 + -1.0 + -1.0 + -1.0)); // -8.0
  EXPECT_EQ(matrix2[0][0][0], -1.0);
  EXPECT_EQ(matrix2[1][0][0], -1.0);
  EXPECT_EQ(matrix2[2][0][0], -1.0);
  EXPECT_EQ(matrix2[0][1][1], -1.0);
  EXPECT_EQ(matrix2[1][1][1], -1.0);
  EXPECT_EQ(matrix2[2][1][1], -1.0);
  EXPECT_EQ(matrix2[0][2][2], -1.0);
  EXPECT_EQ(matrix2[1][2][2], -1.0);
  EXPECT_EQ(matrix2[2][2][2], -1.0);
  EXPECT_EQ(matrix2[0][3][3], -1.0);
  EXPECT_EQ(matrix2[1][3][3], -1.0);
  EXPECT_EQ(matrix2[2][3][3], -1.0);

  return matrix;
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
std::tuple<MatrixPolicy<double, OrderingPolicy>, MatrixPolicy<double, OrderingPolicy>> testMultiMatrixArrayFunction()
{
  // MatrixA: 3x3 with 2 non-zero elements per block
  auto builderA = MatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .WithElement(1, 2)
                      .SetNumberOfBlocks(3);
  // 0 X 0
  // 0 0 X
  // 0 0 0

  // MatrixB: 3x3 with 3 non-zero elements per block
  auto builderB = MatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 0)
                      .WithElement(1, 1)
                      .WithElement(2, 2)
                      .SetNumberOfBlocks(3);
  // X 0 0
  // 0 X 0
  // 0 0 X

  MatrixPolicy<double, OrderingPolicy> matrixA{ builderA };
  MatrixPolicy<double, OrderingPolicy> matrixB{ builderB };

  // Set initial values that differ by blocks
  for (int block = 0; block < 3; ++block)
  {
    matrixA[block][0][1] = static_cast<double>(block + 10);
    matrixA[block][1][2] = static_cast<double>(block * 2 + 20);
    
    matrixB[block][0][0] = static_cast<double>(block * 4);
    matrixB[block][1][1] = static_cast<double>(block * 3);
    matrixB[block][2][2] = static_cast<double>(block * 5);
  }

  // Initial MatrixA values:
  // Block 0: (0,1)=10, (1,2)=20
  // Block 1: (0,1)=11, (1,2)=22
  // Block 2: (0,1)=12, (1,2)=24

  // Initial MatrixB values:
  // Block 0: (0,0)=0, (1,1)=0, (2,2)=0
  // Block 1: (0,0)=4, (1,1)=3, (2,2)=5
  // Block 2: (0,0)=8, (1,1)=6, (2,2)=10

  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& mA, auto&& mB)
    {
      // Use an array function to set element (1,2) in matrixA = element (0,1) in matrixA + element (2,2) in matrixB
      auto tmp = mA.GetBlockVariable();
      mA.ForEachBlock([&](const double& a, const double& b, double& t)
        { t = a + b; },
        mA.GetConstBlockView(0, 1),
        mB.GetConstBlockView(2, 2),
        tmp);
      mA.ForEachBlock([&](const double& t, double& c)
        { c = t; },
        tmp,
        mA.GetBlockView(1, 2));
    }, matrixA, matrixB);

  func(matrixA, matrixB);

  // Check results
  EXPECT_EQ(matrixA[0][1][2], 10 + 0);   // 10
  EXPECT_EQ(matrixA[1][1][2], 11 + 5);   // 16
  EXPECT_EQ(matrixA[2][1][2], 12 + 10);  // 22
  EXPECT_EQ(matrixA[0][0][1], 10.0);
  EXPECT_EQ(matrixA[1][0][1], 11.0);
  EXPECT_EQ(matrixA[2][0][1], 12.0);
  EXPECT_EQ(matrixB[0][0][0], 0.0);
  EXPECT_EQ(matrixB[1][0][0], 4.0);
  EXPECT_EQ(matrixB[2][0][0], 8.0);

  return { matrixA, matrixB };
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
void testMismatchedBlockDimensions()
{
  auto builderA = MatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .WithElement(1, 1)
                      .SetNumberOfBlocks(3);

  auto builderB = MatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .WithElement(1, 1)
                      .SetNumberOfBlocks(4);  // Different number of blocks!

  MatrixPolicy<double, OrderingPolicy> matrixA{ builderA };
  MatrixPolicy<double, OrderingPolicy> matrixB{ builderB };

  // Should succeed at creation (different block counts allowed at creation)
  using SparseMatrixType = MatrixPolicy<double, OrderingPolicy>;
  auto func = SparseMatrixType::Function(
    [](auto&& mA, auto&& mB)
    {
      // This should work when matrices have same block counts
      mA.ForEachBlock([&](const double& a, const double& b, double& c)
        { c = a + b; },
        mA.GetConstBlockView(0, 1),
        mB.GetConstBlockView(0, 1),
        mA.GetBlockView(1, 1));
    }, matrixA, matrixA);
  
  // Should work with matching block counts
  EXPECT_NO_THROW(func(matrixA, matrixA));
  
  // Should throw at invocation when matrices have different block counts
  EXPECT_ANY_THROW(func(matrixA, matrixB));
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
void testMismatchedElementDimensions()
{
  auto builder = MatrixPolicy<double, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(1, 1)
                     .WithElement(2, 2)
                     .SetNumberOfBlocks(3);
  // 0 X 0 0
  // 0 X 0 0
  // 0 0 X 0
  // 0 0 0 0

  MatrixPolicy<double, OrderingPolicy> matrix{ builder };

  // Create the function - this should succeed
  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m)
    {
      // Try to access a block element that doesn't exist
      m.ForEachBlock([&](const double& a, double& b)
        { b = a * 2.0; },
        m.GetConstBlockView(0, 1),
        m.GetBlockView(3, 3));  // Element (3,3) doesn't exist in this sparse matrix
    }, matrix);

  // Should throw when invoking the function because element (3,3) doesn't exist
  EXPECT_ANY_THROW(func(matrix));
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
void testWrongMatrixDimensions()
{
  auto builder1 = MatrixPolicy<double, OrderingPolicy>::Create(4)
                      .WithElement(0, 1)
                      .WithElement(1, 1)
                      .WithElement(2, 2)
                      .WithElement(3, 3)
                      .SetNumberOfBlocks(3);

  auto builder2 = MatrixPolicy<double, OrderingPolicy>::Create(5)  // Different block size (5x5 vs 4x4)
                      .WithElement(0, 1)
                      .WithElement(1, 1)
                      .WithElement(2, 2)
                      .WithElement(3, 3)
                      .SetNumberOfBlocks(3);

  MatrixPolicy<double, OrderingPolicy> matrix1{ builder1 };
  MatrixPolicy<double, OrderingPolicy> matrix2{ builder2 };

  // Create a function with 4x4 block matrix
  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m)
    {
      m.ForEachBlock([&](const double& a, double& b)
        { b = a * 2.0; },
        m.GetConstBlockView(0, 1),
        m.GetBlockView(3, 3));  // Element (3,3) exists in both matrices
    }, matrix1);

  // Should work fine with matrix1
  EXPECT_NO_THROW(func(matrix1));

  // Should also work with matrix2 since it has the same number of blocks (3)
  // and the accessed elements (0,1) and (3,3) exist in its sparsity pattern
  EXPECT_NO_THROW(func(matrix2));
  
  // But if we try to use a matrix with different number of blocks, it should fail
  auto builder3 = MatrixPolicy<double, OrderingPolicy>::Create(4)
                      .WithElement(0, 1)
                      .WithElement(1, 1)
                      .SetNumberOfBlocks(5);  // Different number of blocks!
  MatrixPolicy<double, OrderingPolicy> matrix3{ builder3 };
  
  // Should throw because number of blocks doesn't match (5 vs 3)
  EXPECT_ANY_THROW(func(matrix3));
}

/// @brief Test: Multiple sparse matrices with DIFFERENT block counts from creation (should work)
template<template<class, class> class MatrixPolicy, class OrderingPolicy>
std::tuple<MatrixPolicy<double, OrderingPolicy>, MatrixPolicy<double, OrderingPolicy>> testMultipleSparseMatricesDifferentBlocksFromCreation()
{
  // Create function with matrices having 3 blocks
  auto builder3blocks = MatrixPolicy<double, OrderingPolicy>::Create(4)
                            .WithElement(0, 1)
                            .WithElement(1, 2)
                            .WithElement(2, 3)
                            .SetNumberOfBlocks(3);
  
  MatrixPolicy<double, OrderingPolicy> matrixA_3blocks{ builder3blocks };
  MatrixPolicy<double, OrderingPolicy> matrixB_3blocks{ builder3blocks };
  
  // Initialize 3-block matrices
  for (int block = 0; block < 3; ++block)
  {
    matrixA_3blocks[block][0][1] = static_cast<double>(block + 1);
    matrixA_3blocks[block][1][2] = static_cast<double>(block + 10);
    matrixB_3blocks[block][2][3] = static_cast<double>(block + 100);
  }
  
  // Create function with 3-block matrices
  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& mA, auto&& mB)
    {
      // Compute mA(2,3) = mA(0,1) + mA(1,2) + mB(2,3)
      auto tmp = mA.GetBlockVariable();
      mA.ForEachBlock([&](const double& a, const double& b, const double& c, double& t)
        { t = a + b + c; },
        mA.GetConstBlockView(0, 1),
        mA.GetConstBlockView(1, 2),
        mB.GetConstBlockView(2, 3),
        tmp);
      mA.ForEachBlock([&](const double& t, double& result)
        { result = t; },
        tmp,
        mA.GetBlockView(2, 3));
    }, matrixA_3blocks, matrixB_3blocks);
  
  // Now use with matrices having 4 blocks (different from creation!)
  auto builder4blocks = MatrixPolicy<double, OrderingPolicy>::Create(4)
                            .WithElement(0, 1)
                            .WithElement(1, 2)
                            .WithElement(2, 3)
                            .SetNumberOfBlocks(4);
  
  MatrixPolicy<double, OrderingPolicy> matrixA_4blocks{ builder4blocks };
  MatrixPolicy<double, OrderingPolicy> matrixB_4blocks{ builder4blocks };
  
  // Initialize 4-block matrices
  for (int block = 0; block < 4; ++block)
  {
    matrixA_4blocks[block][0][1] = static_cast<double>(block + 1);
    matrixA_4blocks[block][1][2] = static_cast<double>(block + 10);
    matrixB_4blocks[block][2][3] = static_cast<double>(block + 100);
  }
  
  // Should work with different block count
  EXPECT_NO_THROW(func(matrixA_4blocks, matrixB_4blocks));
  
  // Verify results for first 4 blocks
  EXPECT_EQ(matrixA_4blocks[0][2][3], 1 + 10 + 100);   // 111
  EXPECT_EQ(matrixA_4blocks[1][2][3], 2 + 11 + 101);   // 114
  EXPECT_EQ(matrixA_4blocks[2][2][3], 3 + 12 + 102);   // 117
  EXPECT_EQ(matrixA_4blocks[3][2][3], 4 + 13 + 103);   // 120
  
  return { matrixA_4blocks, matrixB_4blocks };
}

/// @brief Test: Sparse matrix + vector with DIFFERENT block/size from creation (should work)
template<template<class, class> class MatrixPolicy, class OrderingPolicy>
std::tuple<MatrixPolicy<double, OrderingPolicy>, std::vector<double>> testSparseMatrixVectorDifferentBlocksFromCreation()
{
  // Create function with 3-block matrix and 3-element vector
  auto builder3 = MatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .WithElement(1, 2)
                      .SetNumberOfBlocks(3);
  
  MatrixPolicy<double, OrderingPolicy> matrix3{ builder3 };
  std::vector<double> vec3 = { 1.0, 2.0, 3.0 };
  
  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m, auto&& v)
    {
      // m(1,2) = m(0,1) + v
      auto tmp = m.GetBlockVariable();
      m.ForEachBlock([&](const double& a, const double& b, double& t)
        { t = a + b; },
        m.GetConstBlockView(0, 1),
        v,
        tmp);
      m.ForEachBlock([&](const double& t, double& result)
        { result = t; },
        tmp,
        m.GetBlockView(1, 2));
    }, matrix3, vec3);
  
  // Now use with 5-block matrix and 5-element vector (different from creation!)
  auto builder5 = MatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .WithElement(1, 2)
                      .SetNumberOfBlocks(5);
  
  MatrixPolicy<double, OrderingPolicy> matrix5{ builder5 };
  std::vector<double> vec5 = { 10.0, 20.0, 30.0, 40.0, 50.0 };
  
  // Initialize
  for (int block = 0; block < 5; ++block)
  {
    matrix5[block][0][1] = static_cast<double>(block + 1);
  }
  
  // Should work with different block/vector size
  EXPECT_NO_THROW(func(matrix5, vec5));
  
  // Verify results
  EXPECT_EQ(matrix5[0][1][2], 1 + 10);   // 11
  EXPECT_EQ(matrix5[1][1][2], 2 + 20);   // 22
  EXPECT_EQ(matrix5[2][1][2], 3 + 30);   // 33
  EXPECT_EQ(matrix5[3][1][2], 4 + 40);   // 44
  EXPECT_EQ(matrix5[4][1][2], 5 + 50);   // 55
  
  return { matrix5, vec5 };
}

/// @brief Test: MISMATCHED block counts at invocation (should throw)
template<template<class, class> class MatrixPolicy, class OrderingPolicy>
void testMismatchedBlocksAtInvocation()
{
  auto builder3 = MatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .WithElement(1, 2)
                      .SetNumberOfBlocks(3);
  
  auto builder4 = MatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .WithElement(1, 2)
                      .SetNumberOfBlocks(4);
  
  MatrixPolicy<double, OrderingPolicy> matrix3{ builder3 };
  MatrixPolicy<double, OrderingPolicy> matrix4{ builder4 };
  
  // Create function
  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& mA, auto&& mB)
    {
      mA.ForEachBlock([&](const double& a, double& b)
        { b = a * 2.0; },
        mB.GetConstBlockView(0, 1),
        mA.GetBlockView(1, 2));
    }, matrix3, matrix3);
  
  // Should work with matching blocks
  EXPECT_NO_THROW(func(matrix3, matrix3));
  
  // Should throw with mismatched block counts
  EXPECT_ANY_THROW(func(matrix3, matrix4));
  EXPECT_ANY_THROW(func(matrix4, matrix3));
}

/// @brief Test: Multiple sparse matrices with MISMATCHED blocks at invocation (should throw)
template<template<class, class> class MatrixPolicy, class OrderingPolicy>
void testMultipleSparseMatricesMismatchedBlocksAtInvocation()
{
  auto builder3 = MatrixPolicy<double, OrderingPolicy>::Create(4)
                      .WithElement(0, 1)
                      .WithElement(1, 2)
                      .WithElement(2, 3)
                      .SetNumberOfBlocks(3);
  
  auto builder4 = MatrixPolicy<double, OrderingPolicy>::Create(4)
                      .WithElement(0, 1)
                      .WithElement(1, 2)
                      .WithElement(2, 3)
                      .SetNumberOfBlocks(4);
  
  MatrixPolicy<double, OrderingPolicy> matrixA_3{ builder3 };
  MatrixPolicy<double, OrderingPolicy> matrixB_3{ builder3 };
  MatrixPolicy<double, OrderingPolicy> matrixC_4{ builder4 };
  
  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& mA, auto&& mB, auto&& mC)
    {
      mA.ForEachBlock([&](const double& a, const double& b, const double& c, double& result)
        { result = a + b + c; },
        mA.GetConstBlockView(0, 1),
        mB.GetConstBlockView(1, 2),
        mC.GetConstBlockView(2, 3),
        mA.GetBlockView(2, 3));
    }, matrixA_3, matrixB_3, matrixB_3);
  
  // Should work when all have matching blocks
  EXPECT_NO_THROW(func(matrixA_3, matrixB_3, matrixB_3));
  
  // Should throw when one matrix has different block count
  EXPECT_ANY_THROW(func(matrixA_3, matrixB_3, matrixC_4));
  EXPECT_ANY_THROW(func(matrixC_4, matrixB_3, matrixB_3));
}

/// @brief Test: Wrong element structure fails, different blocks succeeds
template<template<class, class> class MatrixPolicy, class OrderingPolicy>
void testWrongStructureAtInvocation()
{
  // Same structure, different blocks
  auto builder3 = MatrixPolicy<double, OrderingPolicy>::Create(4)
                      .WithElement(0, 1)
                      .WithElement(1, 2)
                      .WithElement(2, 3)
                      .SetNumberOfBlocks(3);
  
  auto builder5_same = MatrixPolicy<double, OrderingPolicy>::Create(4)
                           .WithElement(0, 1)
                           .WithElement(1, 2)
                           .WithElement(2, 3)
                           .SetNumberOfBlocks(5);
  
  // Different structure (different elements)
  auto builder5_diff = MatrixPolicy<double, OrderingPolicy>::Create(4)
                           .WithElement(0, 1)
                           .WithElement(1, 1)  // Different!
                           .WithElement(2, 2)  // Different!
                           .SetNumberOfBlocks(5);
  
  MatrixPolicy<double, OrderingPolicy> matrix3{ builder3 };
  MatrixPolicy<double, OrderingPolicy> matrix5_same{ builder5_same };
  MatrixPolicy<double, OrderingPolicy> matrix5_diff{ builder5_diff };
  
  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m)
    {
      m.ForEachBlock([&](const double& a, double& b)
        { b = a * 2.0; },
        m.GetConstBlockView(0, 1),
        m.GetBlockView(1, 2));
    }, matrix3);
  
  // Should work with different block count but same structure
  EXPECT_NO_THROW(func(matrix5_same));
  
  // Should throw with different element structure
  EXPECT_ANY_THROW(func(matrix5_diff));
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<double, OrderingPolicy> testMultipleTemporaries()
{
  auto builder = MatrixPolicy<double, OrderingPolicy>::Create(5)
                     .WithElement(0, 1)
                     .WithElement(1, 2)
                     .WithElement(2, 3)
                     .WithElement(3, 4)
                     .WithElement(4, 4)
                     .SetNumberOfBlocks(4);
  // 0 X 0 0 0
  // 0 0 X 0 0
  // 0 0 0 X 0
  // 0 0 0 0 X
  // 0 0 0 0 X

  MatrixPolicy<double, OrderingPolicy> matrix{ builder };

  // Initialize first two block elements
  for (std::size_t block = 0; block < 4; ++block)
  {
    matrix[block][0][1] = static_cast<double>(block + 1);
    matrix[block][1][2] = static_cast<double>((block + 1) * 10);
  }

  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m)
    {
      // Use TWO temporaries for intermediate calculations
      auto tmp1 = m.GetBlockVariable();
      auto tmp2 = m.GetBlockVariable();

      // tmp1 = (0,1) * (1,2)
      m.ForEachBlock([&](const double& a, const double& b, double& t)
        { t = a * b; },
        m.GetConstBlockView(0, 1),
        m.GetConstBlockView(1, 2),
        tmp1);

      // tmp2 = (0,1) + (1,2)
      m.ForEachBlock([&](const double& a, const double& b, double& t)
        { t = a + b; },
        m.GetConstBlockView(0, 1),
        m.GetConstBlockView(1, 2),
        tmp2);

      // (2,3) = tmp1 + tmp2 (product + sum)
      m.ForEachBlock([&](const double& t1, const double& t2, double& c)
        { c = t1 + t2; },
        tmp1,
        tmp2,
        m.GetBlockView(2, 3));

      // (3,4) = tmp1 - tmp2 (product - sum)
      m.ForEachBlock([&](const double& t1, const double& t2, double& c)
        { c = t1 - t2; },
        tmp1,
        tmp2,
        m.GetBlockView(3, 4));

      // (4,4) = tmp1 * tmp2
      m.ForEachBlock([&](const double& t1, const double& t2, double& c)
        { c = t1 * t2; },
        tmp1,
        tmp2,
        m.GetBlockView(4, 4));
    }, matrix);

  func(matrix);

  // Verify results
  // Block 0: (0,1)=1, (1,2)=10, product=10, sum=11
  EXPECT_EQ(matrix[0][2][3], 10.0 + 11.0);   // 21
  EXPECT_EQ(matrix[0][3][4], 10.0 - 11.0);   // -1
  EXPECT_EQ(matrix[0][4][4], 10.0 * 11.0);   // 110

  // Block 1: (0,1)=2, (1,2)=20, product=40, sum=22
  EXPECT_EQ(matrix[1][2][3], 40.0 + 22.0);   // 62
  EXPECT_EQ(matrix[1][3][4], 40.0 - 22.0);   // 18
  EXPECT_EQ(matrix[1][4][4], 40.0 * 22.0);   // 880

  // Block 3: (0,1)=4, (1,2)=40, product=160, sum=44
  EXPECT_EQ(matrix[3][2][3], 160.0 + 44.0);  // 204
  EXPECT_EQ(matrix[3][3][4], 160.0 - 44.0);  // 116
  EXPECT_EQ(matrix[3][4][4], 160.0 * 44.0);  // 7040

  return matrix;
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<double, OrderingPolicy> testBlockViewReuse()
{
  auto builder = MatrixPolicy<double, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(1, 2)
                     .WithElement(2, 3)
                     .WithElement(3, 3)
                     .SetNumberOfBlocks(3);
  // 0 X 0 0
  // 0 0 X 0
  // 0 0 0 X
  // 0 0 0 X

  MatrixPolicy<double, OrderingPolicy> matrix{ builder };

  for (std::size_t block = 0; block < 3; ++block)
    matrix[block][0][1] = static_cast<double>(block + 1);

  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m)
    {
      // Create block views once
      auto elem01 = m.GetConstBlockView(0, 1);
      auto elem12 = m.GetBlockView(1, 2);
      auto elem23 = m.GetBlockView(2, 3);
      auto elem33 = m.GetBlockView(3, 3);

      // Reuse the same block views in multiple ForEachBlock calls
      // (1,2) = (0,1) * 2
      m.ForEachBlock([&](const double& a, double& b)
        { b = a * 2.0; },
        elem01, elem12);

      // (2,3) = (0,1) + (1,2) (reusing elem01 and elem12)
      m.ForEachBlock([&](const double& a, const double& b, double& c)
        { c = a + b; },
        elem01, elem12, elem23);

      // (3,3) = (2,3) * (1,2) (reusing elem12 and elem23)
      m.ForEachBlock([&](const double& a, const double& b, double& c)
        { c = a * b; },
        elem23, elem12, elem33);
    }, matrix);

  func(matrix);

  // Block 0: (0,1)=1, (1,2)=2, (2,3)=3, (3,3)=6
  EXPECT_EQ(matrix[0][1][2], 2.0);
  EXPECT_EQ(matrix[0][2][3], 3.0);
  EXPECT_EQ(matrix[0][3][3], 6.0);

  // Block 1: (0,1)=2, (1,2)=4, (2,3)=6, (3,3)=24
  EXPECT_EQ(matrix[1][1][2], 4.0);
  EXPECT_EQ(matrix[1][2][3], 6.0);
  EXPECT_EQ(matrix[1][3][3], 24.0);

  // Block 2: (0,1)=3, (1,2)=6, (2,3)=9, (3,3)=54
  EXPECT_EQ(matrix[2][1][2], 6.0);
  EXPECT_EQ(matrix[2][2][3], 9.0);
  EXPECT_EQ(matrix[2][3][3], 54.0);

  return matrix;
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<double, OrderingPolicy> testFunctionReusability()
{
  // Create a function once
  auto builder = MatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .WithElement(2, 2)
                     .SetNumberOfBlocks(2);
  // X 0 0
  // 0 X 0
  // 0 0 X

  MatrixPolicy<double, OrderingPolicy> matrix1{ builder };
  
  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m)
    {
      auto tmp = m.GetBlockVariable();
      m.ForEachBlock([&](const double& a, const double& b, const double& c, double& t)
        { t = a + b + c; },
        m.GetConstBlockView(0, 0),
        m.GetConstBlockView(1, 1),
        m.GetConstBlockView(2, 2),
        tmp);
      m.ForEachBlock([&](double& c, const double& t)
        { c = 2.0 * t; },
        m.GetBlockView(2, 2),
        tmp);
    }, matrix1);

  // Apply to first matrix
  for (std::size_t block = 0; block < 2; ++block)
  {
    matrix1[block][0][0] = static_cast<double>(block);
    matrix1[block][1][1] = static_cast<double>(block + 1);
    matrix1[block][2][2] = static_cast<double>(block + 2);
  }

  func(matrix1);
  EXPECT_EQ(matrix1[0][2][2], 2.0 * (0 + 1 + 2));  // 6
  EXPECT_EQ(matrix1[1][2][2], 2.0 * (1 + 2 + 3));  // 12

  // Apply to second matrix with same dimensions
  MatrixPolicy<double, OrderingPolicy> matrix2{ builder };
  matrix2 = 5.0;
  func(matrix2);
  EXPECT_EQ(matrix2[0][2][2], 2.0 * (5 + 5 + 5));  // 30
  EXPECT_EQ(matrix2[1][2][2], 2.0 * (5 + 5 + 5));  // 30

  // Apply to third matrix with different values
  MatrixPolicy<double, OrderingPolicy> matrix3{ builder };
  matrix3 = 0.0;
  for (std::size_t block = 0; block < 2; ++block)
    matrix3[block][0][0] = static_cast<double>(block * 10);

  func(matrix3);
  EXPECT_EQ(matrix3[0][2][2], 2.0 * (0 + 0 + 0));   // 0
  EXPECT_EQ(matrix3[1][2][2], 2.0 * (10 + 0 + 0));  // 20

  return matrix1;
}

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
std::tuple<MatrixPolicy<double, OrderingPolicy>, MatrixPolicy<double, OrderingPolicy>> testTwoSparseMatricesDifferentStructure()
{
  // MatrixA: 3x3 with specific sparsity pattern, 4 blocks
  auto builderA = MatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .WithElement(1, 1)
                      .WithElement(2, 2)
                      .SetNumberOfBlocks(4);
  // 0 X 0
  // 0 X 0
  // 0 0 X

  // MatrixB: 5x5 with different sparsity pattern and dimensions, but same number of blocks (4)
  auto builderB = MatrixPolicy<double, OrderingPolicy>::Create(5)
                      .WithElement(0, 0)
                      .WithElement(1, 2)
                      .WithElement(2, 3)
                      .WithElement(3, 4)
                      .SetNumberOfBlocks(4);
  // X 0 0 0 0
  // 0 0 X 0 0
  // 0 0 0 X 0
  // 0 0 0 0 X
  // 0 0 0 0 0

  MatrixPolicy<double, OrderingPolicy> matrixA{ builderA };
  MatrixPolicy<double, OrderingPolicy> matrixB{ builderB };

  // Set initial values that differ by blocks
  for (int block = 0; block < 4; ++block)
  {
    matrixA[block][0][1] = static_cast<double>(block * 2 + 1);
    matrixA[block][1][1] = static_cast<double>(block * 3 + 2);
    matrixA[block][2][2] = static_cast<double>(block * 5 + 3);
    
    matrixB[block][0][0] = static_cast<double>(block + 10);
    matrixB[block][1][2] = static_cast<double>(block + 20);
    matrixB[block][2][3] = static_cast<double>(block + 30);
    matrixB[block][3][4] = static_cast<double>(block + 40);
  }

  // Initial MatrixA values by block:
  // Block 0: (0,1)=1, (1,1)=2, (2,2)=3
  // Block 1: (0,1)=3, (1,1)=5, (2,2)=8
  // Block 2: (0,1)=5, (1,1)=8, (2,2)=13
  // Block 3: (0,1)=7, (1,1)=11, (2,2)=18

  // Initial MatrixB values by block:
  // Block 0: (0,0)=10, (1,2)=20, (2,3)=30, (3,4)=40
  // Block 1: (0,0)=11, (1,2)=21, (2,3)=31, (3,4)=41
  // Block 2: (0,0)=12, (1,2)=22, (2,3)=32, (3,4)=42
  // Block 3: (0,0)=13, (1,2)=23, (2,3)=33, (3,4)=43

  auto func = MatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& mA, auto&& mB)
    {
      // Combine data from both matrices despite different structures
      // Set (2,2) in matrixA = (0,1) + (1,1) from matrixA + (0,0) from matrixB
      auto tmp = mA.GetBlockVariable();
      mA.ForEachBlock([&](const double& a, const double& b, const double& c, double& t)
        { t = a + b + c; },
        mA.GetConstBlockView(0, 1),
        mA.GetConstBlockView(1, 1),
        mB.GetConstBlockView(0, 0),
        tmp);
      mA.ForEachBlock([&](const double& t, double& result)
        { result = t; },
        tmp,
        mA.GetBlockView(2, 2));
    }, matrixA, matrixB);

  func(matrixA, matrixB);

  // Check results
  EXPECT_EQ(matrixA[0][2][2], 1 + 2 + 10);    // 13
  EXPECT_EQ(matrixA[1][2][2], 3 + 5 + 11);    // 19
  EXPECT_EQ(matrixA[2][2][2], 5 + 8 + 12);    // 25
  EXPECT_EQ(matrixA[3][2][2], 7 + 11 + 13);   // 31
  
  // Verify other elements unchanged in matrixA
  EXPECT_EQ(matrixA[0][0][1], 1.0);
  EXPECT_EQ(matrixA[1][0][1], 3.0);
  EXPECT_EQ(matrixA[2][0][1], 5.0);
  EXPECT_EQ(matrixA[3][0][1], 7.0);
  
  // Verify matrixB unchanged
  EXPECT_EQ(matrixB[0][0][0], 10.0);
  EXPECT_EQ(matrixB[1][1][2], 21.0);
  EXPECT_EQ(matrixB[2][2][3], 32.0);
  EXPECT_EQ(matrixB[3][3][4], 43.0);

  return { matrixA, matrixB };
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy, template<class> class DenseMatrixPolicy>
std::tuple<SparseMatrixPolicy<double, OrderingPolicy>, DenseMatrixPolicy<double>> testSparseAndDenseMatrixFunction()
{
  // Sparse matrix: 4x4 with some sparsity pattern, 3 blocks
  auto sparseBuilder = SparseMatrixPolicy<double, OrderingPolicy>::Create(4)
                           .WithElement(0, 1)
                           .WithElement(1, 2)
                           .WithElement(2, 3)
                           .WithElement(3, 3)
                           .SetNumberOfBlocks(3);
  // 0 X 0 0
  // 0 0 X 0
  // 0 0 0 X
  // 0 0 0 X

  // Dense matrix: 3 rows (matching number of blocks in sparse matrix) x 4 columns
  // Note: In the implementation, we'll need to validate that the orderings are compatible
  // and that dense.NumRows() == sparse.NumberOfBlocks()

  SparseMatrixPolicy<double, OrderingPolicy> sparseMatrix{ sparseBuilder };
  DenseMatrixPolicy<double> denseMatrix{ 3, 4, 0.0 };

  // Initialize sparse matrix elements
  for (int block = 0; block < 3; ++block)
  {
    sparseMatrix[block][0][1] = static_cast<double>(block + 1);
    sparseMatrix[block][1][2] = static_cast<double>(block + 10);
    sparseMatrix[block][2][3] = static_cast<double>(block + 20);
    sparseMatrix[block][3][3] = static_cast<double>(block + 30);
  }

  // Initialize dense matrix - each row corresponds to a block in the sparse matrix
  for (int row = 0; row < 3; ++row)
  {
    for (int col = 0; col < 4; ++col)
    {
      denseMatrix[row][col] = static_cast<double>(row * 100 + col * 10);
    }
  }

  // Dense matrix values:
  // Row 0: 0, 10, 20, 30
  // Row 1: 100, 110, 120, 130
  // Row 2: 200, 210, 220, 230

  // Create a function that combines sparse blocks with dense rows
  // The function should work with blocks in sparse matrix corresponding to rows in dense matrix
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& sparse, auto&& dense)
    {
      // For each block/row, compute: sparse(3,3) = sparse(0,1) + sparse(1,2) + dense[col0] + dense[col1]
      auto tmp = sparse.GetBlockVariable();
      
      // tmp = sparse(0,1) + sparse(1,2) + dense[0] + dense[1]
      sparse.ForEachBlock([&](const double& s1, const double& s2, const double& d0, const double& d1, double& t)
        { t = s1 + s2 + d0 + d1; },
        sparse.GetConstBlockView(0, 1),
        sparse.GetConstBlockView(1, 2),
        dense.GetConstColumnView(0),  // Dense columns act like sparse block elements
        dense.GetConstColumnView(1),
        tmp);
      
      // sparse(3,3) = tmp
      sparse.ForEachBlock([&](const double& t, double& result)
        { result = t; },
        tmp,
        sparse.GetBlockView(3, 3));
    }, sparseMatrix, denseMatrix);

  func(sparseMatrix, denseMatrix);

  // Check results
  // Block/Row 0: sparse(0,1)=1, sparse(1,2)=10, dense[0]=0, dense[1]=10
  EXPECT_EQ(sparseMatrix[0][3][3], 1 + 10 + 0 + 10);     // 21
  
  // Block/Row 1: sparse(0,1)=2, sparse(1,2)=11, dense[0]=100, dense[1]=110
  EXPECT_EQ(sparseMatrix[1][3][3], 2 + 11 + 100 + 110);  // 223
  
  // Block/Row 2: sparse(0,1)=3, sparse(1,2)=12, dense[0]=200, dense[1]=210
  EXPECT_EQ(sparseMatrix[2][3][3], 3 + 12 + 200 + 210);  // 425

  // Verify other elements unchanged
  EXPECT_EQ(sparseMatrix[0][0][1], 1.0);
  EXPECT_EQ(sparseMatrix[1][1][2], 11.0);
  EXPECT_EQ(sparseMatrix[2][2][3], 22.0);
  
  // Verify dense matrix unchanged
  EXPECT_EQ(denseMatrix[0][0], 0.0);
  EXPECT_EQ(denseMatrix[1][2], 120.0);
  EXPECT_EQ(denseMatrix[2][3], 230.0);

  return { sparseMatrix, denseMatrix };
}

///  @brief Test that mixing incompatible orderings throws an error
/// This should fail: vector-ordered sparse (L>1) mixed with standard-ordered dense (L=1)
template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy, template<class> class DenseMatrixPolicy>
void testIncompatibleOrdering()
{
  // Only run this test if L > 1 (vector ordering)
  static constexpr std::size_t L = OrderingPolicy::GroupVectorSize();
  if constexpr (L > 1)
  {
    // Sparse matrix with vector ordering (L > 1)
    auto sparseBuilder = SparseMatrixPolicy<double, OrderingPolicy>::Create(4)
                             .WithElement(0, 1)
                             .WithElement(1, 2)
                             .SetNumberOfBlocks(3);
    SparseMatrixPolicy<double, OrderingPolicy> sparseMatrix{ sparseBuilder };
    
    // Dense matrix with standard ordering (L = 1)
    DenseMatrixPolicy<double> denseMatrix{ 3, 4, 0.0 };
    
    // This should throw during validation (before lambda execution)
    EXPECT_THROW(
      (SparseMatrixPolicy<double, OrderingPolicy>::Function(
        [](auto&&, auto&&) {},
        sparseMatrix, denseMatrix)),
      std::system_error);
  }
}

/// @brief Test valid combinations of sparse and vector matrices with matching L
template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy, std::size_t L>
std::tuple<SparseMatrixPolicy<double, OrderingPolicy>, micm::VectorMatrix<double, L>> testSparseAndVectorMatrixFunction()
{
  // Verify L matches
  static_assert(OrderingPolicy::GroupVectorSize() == L, "L parameter must match OrderingPolicy GroupVectorSize");
  
  // Sparse matrix: 4x4 with some sparsity pattern, 3 blocks
  auto sparseBuilder = SparseMatrixPolicy<double, OrderingPolicy>::Create(4)
                           .WithElement(0, 1)
                           .WithElement(1, 2)
                           .WithElement(2, 3)
                           .WithElement(3, 3)
                           .SetNumberOfBlocks(3);

  SparseMatrixPolicy<double, OrderingPolicy> sparseMatrix{ sparseBuilder };
  micm::VectorMatrix<double, L> vectorMatrix{ 3, 4, 0.0 };

  // Initialize sparse matrix elements
  for (int block = 0; block < 3; ++block)
  {
    sparseMatrix[block][0][1] = static_cast<double>(block + 1);
    sparseMatrix[block][1][2] = static_cast<double>(block + 10);
    sparseMatrix[block][2][3] = static_cast<double>(block + 20);
    sparseMatrix[block][3][3] = static_cast<double>(block + 30);
  }

  // Initialize vector matrix
  for (int row = 0; row < 3; ++row)
  {
    for (int col = 0; col < 4; ++col)
    {
      vectorMatrix[row][col] = static_cast<double>(row * 100 + col * 10);
    }
  }

  // Create a function that combines sparse blocks with vector matrix rows
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& sparse, auto&& vector)
    {
      auto tmp = sparse.GetBlockVariable();
      
      // tmp = sparse(0,1) + sparse(1,2) + vector[0] + vector[1]
      sparse.ForEachBlock([&](const double& s1, const double& s2, const double& v0, const double& v1, double& t)
        { t = s1 + s2 + v0 + v1; },
        sparse.GetConstBlockView(0, 1),
        sparse.GetConstBlockView(1, 2),
        vector.GetConstColumnView(0),
        vector.GetConstColumnView(1),
        tmp);
      
      // sparse(3,3) = tmp
      sparse.ForEachBlock([&](const double& t, double& result)
        { result = t; },
        tmp,
        sparse.GetBlockView(3, 3));
    }, sparseMatrix, vectorMatrix);

  func(sparseMatrix, vectorMatrix);

  // Check results (same as dense matrix test)
  EXPECT_EQ(sparseMatrix[0][3][3], 1 + 10 + 0 + 10);     // 21
  EXPECT_EQ(sparseMatrix[1][3][3], 2 + 11 + 100 + 110);  // 223
  EXPECT_EQ(sparseMatrix[2][3][3], 3 + 12 + 200 + 210);  // 425

  return { sparseMatrix, vectorMatrix };
}

/// @brief Test that mixing vector matrices with different L values throws an error
template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy, std::size_t DifferentL>
void testIncompatibleVectorOrdering()
{
  static constexpr std::size_t SparseL = OrderingPolicy::GroupVectorSize();
  
  // Only run if L values are different
  if constexpr (SparseL != DifferentL)
  {
    // Sparse matrix with OrderingPolicy L
    auto sparseBuilder = SparseMatrixPolicy<double, OrderingPolicy>::Create(4)
                             .WithElement(0, 1)
                             .WithElement(1, 2)
                             .SetNumberOfBlocks(3);
    SparseMatrixPolicy<double, OrderingPolicy> sparseMatrix{ sparseBuilder };
    
    // Vector matrix with different L
    micm::VectorMatrix<double, DifferentL> vectorMatrix{ 3, 4, 0.0 };
    
    // This should throw during validation (before lambda execution)
    EXPECT_THROW(
      (SparseMatrixPolicy<double, OrderingPolicy>::Function(
        [](auto&&, auto&&) {},
        sparseMatrix, vectorMatrix)),
      std::system_error);
  }
}

/// @brief Test that mixing two sparse matrices with different L values throws an error
template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy1, class OrderingPolicy2>
void testIncompatibleSparseOrdering()
{
  static constexpr std::size_t L1 = OrderingPolicy1::GroupVectorSize();
  static constexpr std::size_t L2 = OrderingPolicy2::GroupVectorSize();
  
  // Only run if L values are different
  if constexpr (L1 != L2)
  {
    // First sparse matrix
    auto builder1 = SparseMatrixPolicy<double, OrderingPolicy1>::Create(4)
                        .WithElement(0, 1)
                        .WithElement(1, 2)
                        .SetNumberOfBlocks(3);
    SparseMatrixPolicy<double, OrderingPolicy1> sparseMatrix1{ builder1 };
    
    // Second sparse matrix with different ordering
    auto builder2 = SparseMatrixPolicy<double, OrderingPolicy2>::Create(4)
                        .WithElement(0, 1)
                        .WithElement(1, 2)
                        .SetNumberOfBlocks(3);
    SparseMatrixPolicy<double, OrderingPolicy2> sparseMatrix2{ builder2 };
    
    // This should throw during validation (before lambda execution)
    // Use a simple lambda that doesn't access mismatched data to avoid template instantiation issues
    EXPECT_THROW(
      (SparseMatrixPolicy<double, OrderingPolicy1>::Function(
        [](auto&&, auto&&) {}, 
        sparseMatrix1, sparseMatrix2)),
      std::system_error);
  }
}

/// @brief Test const-correctness with sparse matrices
template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testConstSparseMatrixFunction()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(1, 2)
                     .WithElement(2, 3)
                     .WithElement(3, 3)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  
  // Set initial values
  for (int block = 0; block < 3; ++block)
  {
    matrix[block][0][1] = static_cast<double>(block + 1);
    matrix[block][1][2] = static_cast<double>(block + 10);
  }
  
  // Create a const reference
  const SparseMatrixPolicy<double, OrderingPolicy>& const_matrix = matrix;
  
  // Create a function that only reads from the matrix
  auto read_func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m)
    {
      auto tmp = m.GetBlockVariable();
      // Only use GetConstBlockView - should work with const matrices
      m.ForEachBlock([&](const double& a, const double& b, double& t)
        { t = a + b; },
        m.GetConstBlockView(0, 1),
        m.GetConstBlockView(1, 2),
        tmp);
      
      // Verify we can read the values (no writes to m)
      double sum = 0.0;
      m.ForEachBlock([&sum](const double& val)
        { sum += val; },
        m.GetConstBlockView(2, 3));
    }, const_matrix);
  
  // Should work fine with const matrix
  EXPECT_NO_THROW(read_func(const_matrix));
  
  // Verify original matrix unchanged
  EXPECT_EQ(matrix[0][0][1], 1.0);
  EXPECT_EQ(matrix[1][1][2], 11.0);
}

/// @brief Test edge case with empty sparse matrices
template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testEmptySparseMatrixFunction()
{
  // Test with 0 non-zero elements
  auto empty_builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(4)
                           .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> empty_matrix{ empty_builder };
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m)
    {
      // With 0 non-zero elements, GetBlockView would throw
      // So just iterate (0 blocks will be accessed)
      auto tmp = m.GetBlockVariable();
      // This inner lambda never executes if there are no blocks
    }, empty_matrix);
  
  // Should not throw, function executes but has nothing to do
  EXPECT_NO_THROW(func(empty_matrix));
  
  // Test with 0 blocks
  auto zero_blocks_builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(4)
                                 .WithElement(0, 1)
                                 .SetNumberOfBlocks(0);
  
  SparseMatrixPolicy<double, OrderingPolicy> zero_blocks{ zero_blocks_builder };
  
  auto func2 = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&&)
    {
      // Should never iterate
    }, zero_blocks);
  
  EXPECT_NO_THROW(func2(zero_blocks));
}

// ============================================================================
// Vector Support Tests for Sparse Matrices
// ============================================================================

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
SparseMatrixPolicy<double, OrderingPolicy> testVectorInSparseMatrixFunction()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(1, 2)
                     .WithElement(2, 3)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<double> vec = { 5.0, 10.0, 15.0 };
  
  // Set matrix values
  for (int block = 0; block < 3; ++block)
  {
    matrix[block][0][1] = static_cast<double>(block + 1);
    matrix[block][1][2] = static_cast<double>(block + 2);
    matrix[block][2][3] = static_cast<double>(block + 3);
  }
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m, auto&& v)
    {
      // Read from vector, write to matrix
      m.ForEachBlock([&](const double& in_vec, double& out_mat)
        { out_mat = in_vec * 2.0; },
        v,
        m.GetBlockView(0, 1));
    }, matrix, vec);
  
  func(matrix, vec);
  
  // Check results
  EXPECT_EQ(matrix[0][0][1], 10.0);  // vec[0] * 2
  EXPECT_EQ(matrix[1][0][1], 20.0);  // vec[1] * 2
  EXPECT_EQ(matrix[2][0][1], 30.0);  // vec[2] * 2
  
  return matrix;
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testVectorTooSmall()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<double> vec = { 1.0, 2.0 };  // Too small - needs 3 elements
  
  // Should succeed at creation (block counts can differ at creation)
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
      [](auto&& m, auto&& v) {
        m.ForEachBlock([&](const double& a, double& b)
          { b = a; },
          v,
          m.GetBlockView(0, 1));
      }, matrix, vec);
  
  // Should throw at invocation when vector size doesn't match block count
  EXPECT_ANY_THROW(func(matrix, vec));
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testVectorTooLarge()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<double> vec = { 1.0, 2.0, 3.0, 4.0 };  // Too large
  
  // Should succeed at creation (block counts can differ at creation)
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
      [](auto&& m, auto&& v) {
        m.ForEachBlock([&](const double& a, double& b)
          { b = a; },
          v,
          m.GetBlockView(0, 1));
      }, matrix, vec);
  
  // Should throw at invocation when vector size doesn't match block count
  EXPECT_ANY_THROW(func(matrix, vec));
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testEmptyVectorNonEmptySparseMatrix()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<double> vec;  // Empty
  
  // Should succeed at creation
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
      [](auto&& m, auto&& v) {
        m.ForEachBlock([&](const double& a, double& b)
          { b = a; },
          v,
          m.GetBlockView(0, 1));
      }, matrix, vec);
  
  // Should throw at invocation when vector size doesn't match
  EXPECT_ANY_THROW(func(matrix, vec));
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testNonEmptyVectorEmptySparseMatrix()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(0);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<double> vec = { 1.0, 2.0, 3.0 };
  
  // Should succeed at creation
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
      [](auto&& m, auto&& v) {
        m.ForEachBlock([&](const double& a, double& b)
          { b = a; },
          v,
          m.GetBlockView(0, 1));
      }, matrix, vec);
  
  // Should throw at invocation when vector size doesn't match block count
  EXPECT_ANY_THROW(func(matrix, vec));
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testEmptyVectorEmptySparseMatrix()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(0);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<double> vec;  // Empty
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
      [](auto&&, auto&&) {
        // Should never execute
        FAIL() << "Function should not be called for empty matrix/vector";
      }, matrix, vec);
  
  // Should execute without error (zero iterations)
  EXPECT_NO_THROW(func(matrix, vec));
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testMultipleVectorsDifferentSizes()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<double> vec1 = { 1.0, 2.0, 3.0 };
  std::vector<double> vec2 = { 4.0, 5.0 };  // Different size
  
  // Should succeed at creation (different block counts allowed at creation)
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
      [](auto&& m, auto&& v1, auto&& v2) {
        m.ForEachBlock([&](const double& a, const double& b, double& c)
          { c = a + b; },
          v1,
          v2,
          m.GetBlockView(0, 1));
      }, matrix, vec1, vec2);
  
  // Should throw at invocation because vectors have different sizes
  EXPECT_ANY_THROW(func(matrix, vec1, vec2));
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testMultipleVectorsSameSize()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<double> vec1 = { 1.0, 2.0, 3.0 };
  std::vector<double> vec2 = { 4.0, 5.0, 6.0 };
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m, auto&& v1, auto&& v2)
    {
      m.ForEachBlock([&](const double& a, const double& b, double& out)
        { out = a + b; },
        v1, v2, m.GetBlockView(0, 1));
    }, matrix, vec1, vec2);
  
  func(matrix, vec1, vec2);
  
  EXPECT_EQ(matrix[0][0][1], 5.0);   // 1 + 4
  EXPECT_EQ(matrix[1][0][1], 7.0);   // 2 + 5
  EXPECT_EQ(matrix[2][0][1], 9.0);   // 3 + 6
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testMultipleSparseMatricesOneVector()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix1{ builder };
  SparseMatrixPolicy<double, OrderingPolicy> matrix2{ builder };
  std::vector<double> vec = { 10.0, 20.0, 30.0 };
  
  for (int block = 0; block < 3; ++block)
  {
    matrix1[block][0][1] = static_cast<double>(block + 1);
    matrix2[block][0][1] = 0.0;
  }
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m1, auto&& m2, auto&& v)
    {
      m1.ForEachBlock([&](const double& a, const double& b, double& out)
        { out = a + b; },
        m1.GetConstBlockView(0, 1),
        v,
        m2.GetBlockView(0, 1));
    }, matrix1, matrix2, vec);
  
  func(matrix1, matrix2, vec);
  
  EXPECT_EQ(matrix2[0][0][1], 11.0);  // 1 + 10
  EXPECT_EQ(matrix2[1][0][1], 22.0);  // 2 + 20
  EXPECT_EQ(matrix2[2][0][1], 33.0);  // 3 + 30
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testMultipleSparseMatricesDifferentBlocksVector()
{
  auto builder1 = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .SetNumberOfBlocks(3);
  auto builder2 = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .SetNumberOfBlocks(4);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix1{ builder1 };
  SparseMatrixPolicy<double, OrderingPolicy> matrix2{ builder2 };
  std::vector<double> vec = { 1.0, 2.0, 3.0 };
  
  // Should succeed at creation (different block counts allowed at creation)
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
      [](auto&& m1, auto&& m2, auto&& v) {
        m1.ForEachBlock([&](const double& a, double& b)
          { b = a; },
          v,
          m1.GetBlockView(0, 1));
      }, matrix1, matrix2, vec);
  
  // Should throw at invocation because matrices have different block counts
  EXPECT_ANY_THROW(func(matrix1, matrix2, vec));
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testVectorSizeMatchesOneSparseMatrixOnly()
{
  auto builder1 = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .SetNumberOfBlocks(3);
  auto builder2 = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                      .WithElement(0, 1)
                      .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix1{ builder1 };
  SparseMatrixPolicy<double, OrderingPolicy> matrix2{ builder2 };
  std::vector<double> vec1 = { 1.0, 2.0, 3.0 };
  std::vector<double> vec2 = { 4.0, 5.0 };  // Wrong size
  
  // Should succeed at creation (different block counts allowed at creation)
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
      [](auto&& m1, auto&& m2, auto&& v1, auto&& v2) {
        m1.ForEachBlock([&](const double& a, double& b)
          { b = a; },
          v1,
          m1.GetBlockView(0, 1));
      }, matrix1, matrix2, vec1, vec2);
  
  // Should throw at invocation because vectors have different sizes
  EXPECT_ANY_THROW(func(matrix1, matrix2, vec1, vec2));
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
SparseMatrixPolicy<double, OrderingPolicy> testConstVectorSparse()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  const std::vector<double> vec = { 5.0, 10.0, 15.0 };
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m, auto&& v)
    {
      m.ForEachBlock([&](const double& a, double& b)
        { b = a * 3.0; },
        v,
        m.GetBlockView(0, 1));
    }, matrix, vec);
  
  func(matrix, vec);
  
  EXPECT_EQ(matrix[0][0][1], 15.0);
  EXPECT_EQ(matrix[1][0][1], 30.0);
  EXPECT_EQ(matrix[2][0][1], 45.0);
  
  return matrix;
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
std::tuple<SparseMatrixPolicy<double, OrderingPolicy>, std::vector<double>> testMutableVectorSparse()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<double> vec = { 5.0, 10.0, 15.0 };
  
  for (int block = 0; block < 3; ++block)
    matrix[block][0][1] = static_cast<double>(block + 1);
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m, auto&& v)
    {
      m.ForEachBlock([&](const double& a, double& b)
        { b = a * 3.0; },
        m.GetConstBlockView(0, 1),
        v);
    }, matrix, vec);
  
  func(matrix, vec);
  
  // Vector should be modified
  EXPECT_EQ(vec[0], 3.0);   // 1 * 3
  EXPECT_EQ(vec[1], 6.0);   // 2 * 3
  EXPECT_EQ(vec[2], 9.0);   // 3 * 3
  
  return std::make_tuple(matrix, vec);
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testFunctionReusabilityWithVectorsSparse()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix1{ builder };
  SparseMatrixPolicy<double, OrderingPolicy> matrix2{ builder };
  std::vector<double> vec1 = { 1.0, 2.0, 3.0 };
  std::vector<double> vec2 = { 10.0, 20.0, 30.0 };
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m, auto&& v)
    {
      m.ForEachBlock([&](const double& a, double& b)
        { b = a * 2.0; },
        v,
        m.GetBlockView(0, 1));
    }, matrix1, vec1);
  
  func(matrix1, vec1);
  EXPECT_EQ(matrix1[0][0][1], 2.0);
  EXPECT_EQ(matrix1[1][0][1], 4.0);
  EXPECT_EQ(matrix1[2][0][1], 6.0);
  
  func(matrix2, vec2);
  EXPECT_EQ(matrix2[0][0][1], 20.0);
  EXPECT_EQ(matrix2[1][0][1], 40.0);
  EXPECT_EQ(matrix2[2][0][1], 60.0);
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testFunctionInvocationWithWrongSizedVectorSparse()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<double> vec1 = { 1.0, 2.0, 3.0 };
  std::vector<double> vec2 = { 1.0, 2.0 };  // Wrong size
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
      [](auto&&, auto&&) {}, matrix, vec1);
  
  EXPECT_THROW(
      try { func(matrix, vec2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::InvalidVector));
        throw;
      },
      std::system_error);
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testArraySupportSparse()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::array<double, 3> arr = { 7.0, 14.0, 21.0 };
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m, auto&& a)
    {
      m.ForEachBlock([&](const double& val, double& out)
        { out = val / 7.0; },
        a,
        m.GetBlockView(0, 1));
    }, matrix, arr);
  
  func(matrix, arr);
  
  EXPECT_EQ(matrix[0][0][1], 1.0);
  EXPECT_EQ(matrix[1][0][1], 2.0);
  EXPECT_EQ(matrix[2][0][1], 3.0);
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testMixedVectorBlockViewBlockVariable()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(1, 2)
                     .WithElement(2, 3)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<double> vec = { 100.0, 200.0, 300.0 };
  
  for (int block = 0; block < 3; ++block)
  {
    matrix[block][0][1] = static_cast<double>(block + 1);
    matrix[block][1][2] = static_cast<double>(block + 2);
    matrix[block][2][3] = static_cast<double>(block + 3);
  }
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m, auto&& v)
    {
      auto tmp = m.GetBlockVariable();
      m.ForEachBlock([&](const double& a, const double& b, const double& c, double& t)
        { t = a + b + c; },
        m.GetConstBlockView(0, 1),
        m.GetConstBlockView(1, 2),
        v,
        tmp);
      m.ForEachBlock([&](const double& t, double& out)
        { out = t * 2.0; },
        tmp,
        m.GetBlockView(2, 3));
    }, matrix, vec);
  
  func(matrix, vec);
  
  EXPECT_EQ(matrix[0][2][3], 206.0);  // (1 + 2 + 100) * 2
  EXPECT_EQ(matrix[1][2][3], 410.0);  // (1 + 2 + 200) * 2
  EXPECT_EQ(matrix[2][2][3], 614.0);  // (1 + 2 + 300) * 2
}

template<template<class, class> class SparseMatrixPolicy, class OrderingPolicy>
void testIntegerVectorSparse()
{
  auto builder = SparseMatrixPolicy<double, OrderingPolicy>::Create(3)
                     .WithElement(0, 1)
                     .SetNumberOfBlocks(3);
  
  SparseMatrixPolicy<double, OrderingPolicy> matrix{ builder };
  std::vector<int> int_vec = { 5, 10, 15 };
  
  auto func = SparseMatrixPolicy<double, OrderingPolicy>::Function(
    [](auto&& m, auto&& v)
    {
      m.ForEachBlock([&](const int& i, double& out)
        { out = static_cast<double>(i) * 1.5; },
        v,
        m.GetBlockView(0, 1));
    }, matrix, int_vec);
  
  func(matrix, int_vec);
  
  EXPECT_DOUBLE_EQ(matrix[0][0][1], 7.5);
  EXPECT_DOUBLE_EQ(matrix[1][0][1], 15.0);
  EXPECT_DOUBLE_EQ(matrix[2][0][1], 22.5);
}

