#include <micm/util/sparse_matrix.hpp>

#include <gtest/gtest.h>

template<template<class, class> class MatrixPolicy, class OrderingPolicy>
MatrixPolicy<double, OrderingPolicy> testZeroMatrix()
{
  auto builder = MatrixPolicy<double, OrderingPolicy>::Create(3);
  auto row_ids = builder.RowIdsVector();
  auto row_starts = builder.RowStartVector();

  EXPECT_EQ(builder.NumberOfElements(), 0);
  EXPECT_EQ(row_ids.size(), 0);
  EXPECT_EQ(row_starts.size(), 4);
  EXPECT_EQ(row_starts[0], 0);
  EXPECT_EQ(row_starts[1], 0);
  EXPECT_EQ(row_starts[2], 0);
  EXPECT_EQ(row_starts[3], 0);

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
  auto row_ids = builder.RowIdsVector();
  auto row_starts = builder.RowStartVector();

  EXPECT_EQ(builder.NumberOfElements(), 0);
  EXPECT_EQ(row_ids.size(), 0);
  EXPECT_EQ(row_starts.size(), 4);
  EXPECT_EQ(row_starts[0], 0);
  EXPECT_EQ(row_starts[1], 0);
  EXPECT_EQ(row_starts[2], 0);
  EXPECT_EQ(row_starts[3], 0);

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
MatrixPolicy<int, OrderingPolicy> testSingleBlockMatrix()
{
  auto builder = MatrixPolicy<int, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(3, 2)
                     .WithElement(0, 1)
                     .WithElement(2, 3)
                     .WithElement(2, 1);
  // 0 X 0 0
  // 0 0 0 0
  // 0 X 0 X
  // 0 0 X 0
  {
    auto row_ids = builder.RowIdsVector();
    auto row_starts = builder.RowStartVector();

    EXPECT_EQ(builder.NumberOfElements(), 4);
    EXPECT_EQ(row_ids.size(), 4);
    EXPECT_EQ(row_ids[0], 1);
    EXPECT_EQ(row_ids[1], 1);
    EXPECT_EQ(row_ids[2], 3);
    EXPECT_EQ(row_ids[3], 2);
    EXPECT_EQ(row_starts.size(), 5);
    EXPECT_EQ(row_starts[0], 0);
    EXPECT_EQ(row_starts[1], 1);
    EXPECT_EQ(row_starts[2], 1);
    EXPECT_EQ(row_starts[3], 3);
    EXPECT_EQ(row_starts[4], 4);
  }

  MatrixPolicy<int, OrderingPolicy> matrix{ builder };

  {
    auto& row_ids = matrix.RowIdsVector();
    auto& row_starts = matrix.RowStartVector();
    EXPECT_EQ(row_ids.size(), 4);
    EXPECT_EQ(row_ids[0], 1);
    EXPECT_EQ(row_ids[1], 1);
    EXPECT_EQ(row_ids[2], 3);
    EXPECT_EQ(row_ids[3], 2);
    EXPECT_EQ(row_starts.size(), 5);
    EXPECT_EQ(row_starts[0], 0);
    EXPECT_EQ(row_starts[1], 1);
    EXPECT_EQ(row_starts[2], 1);
    EXPECT_EQ(row_starts[3], 3);
    EXPECT_EQ(row_starts[4], 4);
  }

  EXPECT_EQ(matrix.FlatBlockSize(), 4);

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
      try { std::size_t elem = matrix.VectorIndex(1, 1); } catch (const std::system_error& e) {
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
      try { matrix[0][1][1] = 2; } catch (const std::system_error& e) {
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
                     .WithElement(0, 1)
                     .WithElement(2, 3)
                     .WithElement(2, 1);
  // 0 X 0 0
  // 0 0 0 0
  // 0 X 0 X
  // 0 0 X 0
  MatrixPolicy<int, OrderingPolicy> orig_matrix{ builder };
  orig_matrix[0][2][1] = 45;
  orig_matrix[0][3][2] = 42;
  orig_matrix[0][2][3] = 21;
  auto& row_ids = orig_matrix.RowIdsVector();
  EXPECT_EQ(row_ids.size(), 4);
  const MatrixPolicy<int, OrderingPolicy> matrix = orig_matrix;

  {
    auto& row_ids = matrix.RowIdsVector();
    auto& row_starts = matrix.RowStartVector();
    EXPECT_EQ(row_ids.size(), 4);
    EXPECT_EQ(row_ids[0], 1);
    EXPECT_EQ(row_ids[1], 1);
    EXPECT_EQ(row_ids[2], 3);
    EXPECT_EQ(row_ids[3], 2);
    EXPECT_EQ(row_starts.size(), 5);
    EXPECT_EQ(row_starts[0], 0);
    EXPECT_EQ(row_starts[1], 1);
    EXPECT_EQ(row_starts[2], 1);
    EXPECT_EQ(row_starts[3], 3);
    EXPECT_EQ(row_starts[4], 4);
  }

  EXPECT_EQ(matrix.FlatBlockSize(), 4);

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
      try { std::size_t elem = matrix.VectorIndex(1, 1); } catch (const std::system_error& e) {
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
MatrixPolicy<int, OrderingPolicy> testMultiBlockMatrix()
{
  auto builder = MatrixPolicy<int, OrderingPolicy>::Create(4)
                     .WithElement(0, 1)
                     .WithElement(3, 2)
                     .WithElement(0, 1)
                     .WithElement(2, 3)
                     .WithElement(2, 1)
                     .InitialValue(24)
                     .SetNumberOfBlocks(3);
  // 0 X 0 0
  // 0 0 0 0
  // 0 X 0 X
  // 0 0 X 0
  auto row_ids = builder.RowIdsVector();
  auto row_starts = builder.RowStartVector();

  EXPECT_EQ(builder.NumberOfElements(), 4 * 3);
  EXPECT_EQ(row_ids.size(), 4);
  EXPECT_EQ(row_ids[0], 1);
  EXPECT_EQ(row_ids[1], 1);
  EXPECT_EQ(row_ids[2], 3);
  EXPECT_EQ(row_ids[3], 2);
  EXPECT_EQ(row_starts.size(), 5);
  EXPECT_EQ(row_starts[0], 0);
  EXPECT_EQ(row_starts[1], 1);
  EXPECT_EQ(row_starts[2], 1);
  EXPECT_EQ(row_starts[3], 3);
  EXPECT_EQ(row_starts[4], 4);

  MatrixPolicy<int, OrderingPolicy> matrix{ builder };

  EXPECT_EQ(matrix.FlatBlockSize(), 4);

  EXPECT_EQ(matrix[0][2][1], 24);
  matrix[0][2][1] = 45;
  EXPECT_EQ(matrix[0][2][1], 45);

  matrix[2][2][3] = 64;
  EXPECT_EQ(matrix[2][2][3], 64);

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
      try { std::size_t elem = matrix.VectorIndex(4, 0, 2); } catch (const std::system_error& e) {
        EXPECT_EQ(e.code().value(), static_cast<int>(MicmMatrixErrc::ElementOutOfRange));
        throw;
      },
      std::system_error);
  EXPECT_THROW(
      try { std::size_t elem = matrix.VectorIndex(1, 1, 1); } catch (const std::system_error& e) {
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
      try { matrix[3][0][0] = 2; } catch (const std::system_error& e) {
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
      try { matrix[0][1][1] = 2; } catch (const std::system_error& e) {
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
