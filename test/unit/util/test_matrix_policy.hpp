// Tests of common matrix functions
#include <micm/util/reducers.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestSmallMatrix()
{
  MatrixPolicy<micm::Real> matrix(3, 5);

  matrix[1][3] = 64.7;
  matrix[0][0] = 41.2;
  matrix[2][4] = 102.3;

  EXPECT_EQ(matrix[1][3], static_cast<micm::Real>(64.7));
  EXPECT_EQ(matrix[0][0], static_cast<micm::Real>(41.2));
  EXPECT_EQ(matrix[2][4], static_cast<micm::Real>(102.3));

  std::vector<micm::Real>& data = matrix.AsVector();

  EXPECT_GE(data.size(), 3);

  return matrix;
}

template<template<class> class MatrixPolicy>
const MatrixPolicy<micm::Real> TestSmallConstMatrix()
{
  MatrixPolicy<micm::Real> matrix(3, 5);

  matrix[1][3] = 64.7;
  matrix[0][0] = 41.2;
  matrix[2][4] = 102.3;

  const MatrixPolicy<micm::Real> const_matrix = matrix;

  EXPECT_EQ(const_matrix[1][3], static_cast<micm::Real>(64.7));
  EXPECT_EQ(const_matrix[0][0], static_cast<micm::Real>(41.2));
  EXPECT_EQ(const_matrix[2][4], static_cast<micm::Real>(102.3));

  const std::vector<micm::Real>& data = const_matrix.AsVector();

  EXPECT_GE(data.size(), 3);

  return const_matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestInializeMatrix()
{
  MatrixPolicy<micm::Real> matrix{ 2, 3, 12.4 };

  EXPECT_EQ(matrix[0][0], static_cast<micm::Real>(12.4));
  EXPECT_EQ(matrix[1][0], static_cast<micm::Real>(12.4));
  EXPECT_EQ(matrix[1][2], static_cast<micm::Real>(12.4));

  return matrix;
}

template<template<class> class MatrixPolicy>
const MatrixPolicy<micm::Real> TestInializeConstMatrix()
{
  const MatrixPolicy<micm::Real> matrix{ 2, 3, 12.4 };

  EXPECT_EQ(matrix[0][0], static_cast<micm::Real>(12.4));
  EXPECT_EQ(matrix[1][0], static_cast<micm::Real>(12.4));
  EXPECT_EQ(matrix[1][2], static_cast<micm::Real>(12.4));

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<int> TestLoopOverMatrix()
{
  MatrixPolicy<int> matrix(3, 4, 0);
  for (micm::Index i{}; i < matrix.NumRows(); ++i)
  {
    for (micm::Index j{}; j < matrix.NumColumns(); ++j)
    {
      matrix[i][j] = i * 100 + j;
    }
  }

  EXPECT_EQ(matrix[0][0], 0);
  EXPECT_EQ(matrix[1][2], 102);
  EXPECT_EQ(matrix[2][3], 203);
  EXPECT_EQ(matrix[0][3], 3);

  return matrix;
}

template<template<class> class MatrixPolicy>
const MatrixPolicy<int> TestLoopOverConstMatrix()
{
  MatrixPolicy<int> matrix(3, 4, 0);
  for (micm::Index i{}; i < matrix.NumRows(); ++i)
  {
    for (micm::Index j{}; j < matrix.NumColumns(); ++j)
    {
      matrix[i][j] = i * 100 + j;
    }
  }

  const MatrixPolicy<int> const_matrix = matrix;

  EXPECT_EQ(const_matrix[0][0], 0);
  EXPECT_EQ(const_matrix[1][2], 102);
  EXPECT_EQ(const_matrix[2][3], 203);
  EXPECT_EQ(const_matrix[0][3], 3);

  return const_matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<int> TestStrides()
{
  MatrixPolicy<int> matrix(3, 4, 0);

  for (micm::Index i = 0; i < matrix.NumRows(); ++i)
  {
    for (micm::Index j = 0; j < matrix.NumColumns(); ++j)
    {
      matrix.AsVector()[i * matrix.RowStride() + j * matrix.ColumnStride()] = i * 100 + j;
    }
  }

  EXPECT_EQ(matrix[0][0], 0);
  EXPECT_EQ(matrix[0][1], 1);
  EXPECT_EQ(matrix[0][2], 2);
  EXPECT_EQ(matrix[0][3], 3);
  EXPECT_EQ(matrix[1][0], 100);
  EXPECT_EQ(matrix[1][1], 101);
  EXPECT_EQ(matrix[1][2], 102);
  EXPECT_EQ(matrix[1][3], 103);
  EXPECT_EQ(matrix[2][0], 200);
  EXPECT_EQ(matrix[2][1], 201);
  EXPECT_EQ(matrix[2][2], 202);
  EXPECT_EQ(matrix[2][3], 203);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestConversionToVector()
{
  MatrixPolicy<micm::Real> matrix{ 2, 3, 0.0 };

  matrix[1][0] = 13.2;
  matrix[1][1] = 31.2;
  matrix[1][2] = 314.2;

  std::vector<micm::Real> slice = matrix[1];

  EXPECT_EQ(slice[0], static_cast<micm::Real>(13.2));
  EXPECT_EQ(slice[1], static_cast<micm::Real>(31.2));
  EXPECT_EQ(slice[2], static_cast<micm::Real>(314.2));

  return matrix;
}

template<template<class> class MatrixPolicy>
const MatrixPolicy<micm::Real> TestConstConversionToVector()
{
  MatrixPolicy<micm::Real> matrix{ 2, 3, 0.0 };

  matrix[1][0] = 13.2;
  matrix[1][1] = 31.2;
  matrix[1][2] = 314.2;

  const MatrixPolicy<micm::Real> const_matrix = matrix;
  std::vector<micm::Real> slice = const_matrix[1];

  EXPECT_EQ(slice[0], static_cast<micm::Real>(13.2));
  EXPECT_EQ(slice[1], static_cast<micm::Real>(31.2));
  EXPECT_EQ(slice[2], static_cast<micm::Real>(314.2));

  return const_matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestConversionFromVector()
{
  MatrixPolicy<micm::Real> zero_matrix = std::vector<std::vector<micm::Real>>{};

  EXPECT_EQ(zero_matrix.NumRows(), 0);

  std::vector<std::vector<micm::Real>> vec = { { 412.3, 32.4, 41.3 }, { 5.33, -0.3, 31.2 } };

  MatrixPolicy<micm::Real> matrix = vec;

  EXPECT_EQ(matrix.NumRows(), 2);
  EXPECT_EQ(matrix.NumColumns(), 3);
  EXPECT_EQ(matrix[0].Size(), 3);
  EXPECT_EQ(matrix[0][0], static_cast<micm::Real>(412.3));
  EXPECT_EQ(matrix[0][1], static_cast<micm::Real>(32.4));
  EXPECT_EQ(matrix[0][2], static_cast<micm::Real>(41.3));
  EXPECT_EQ(matrix[1].Size(), 3);
  EXPECT_EQ(matrix[1][0], static_cast<micm::Real>(5.33));
  EXPECT_EQ(matrix[1][1], static_cast<micm::Real>(-0.3));
  EXPECT_EQ(matrix[1][2], static_cast<micm::Real>(31.2));

  std::vector<std::vector<int>> bad_vector = { { 3 }, { 4, 5 }, { 5 } };

  MatrixPolicy<int> bad_matrix;
  EXPECT_ANY_THROW(bad_matrix = bad_vector);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestAssignmentFromVector()
{
  MatrixPolicy<micm::Real> matrix{ 4, 3, 0.0 };
  std::vector<micm::Real> other = { 12.3, 15.1, 24.3 };
  std::vector<micm::Real> big_other = { 14.3, 52.3, 65.7, 16.34 };
  std::vector<micm::Real> small_other = { 13.2, 52.8 };

  matrix[2] = other;

  EXPECT_EQ(matrix[0][0], 0.0);
  EXPECT_EQ(matrix[2][0], static_cast<micm::Real>(12.3));
  EXPECT_EQ(matrix[2][1], static_cast<micm::Real>(15.1));
  EXPECT_EQ(matrix[2][2], static_cast<micm::Real>(24.3));
  EXPECT_EQ(matrix[3][0], 0.0);

  matrix[2] = big_other;

  EXPECT_EQ(matrix[0][0], 0.0);
  EXPECT_EQ(matrix[2][0], static_cast<micm::Real>(14.3));
  EXPECT_EQ(matrix[2][1], static_cast<micm::Real>(52.3));
  EXPECT_EQ(matrix[2][2], static_cast<micm::Real>(65.7));
  EXPECT_EQ(matrix[3][0], 0.0);

  EXPECT_ANY_THROW(matrix[2] = small_other);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestAxpy()
{
  micm::Index num_rows = 4;
  micm::Index num_columns = 3;
  MatrixPolicy<micm::Real> matrix{ num_rows, num_columns, 100.0 };
  MatrixPolicy<micm::Real> other{ num_rows, num_columns, 200.0 };
  micm::Real alpha = 1.39;
  micm::Real sum = 0.0;
  micm::Real result = 0.0;

  for (micm::Index i = 0; i < num_rows; ++i)
  {
    for (micm::Index j = 0; j < num_columns; ++j)
    {
      auto y = i * 10.3 + j * 100.5;
      auto x = i * 1.7 + j * 10.2;
      matrix[i][j] = y;
      other[i][j] = x;
      sum += y + alpha * x;
    }
  }

  matrix.CopyToDevice();
  other.CopyToDevice();
  matrix.Axpy(alpha, other);
  matrix.CopyToHost();

  for (micm::Index i = 0; i < num_rows; ++i)
  {
    for (micm::Index j = 0; j < num_columns; ++j)
    {
      result += matrix[i][j];
    }
  }
  EXPECT_NEAR(sum, result, 1.0e-5);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestForEach()
{
  MatrixPolicy<micm::Real> matrix{ 4, 3, 0.0 };
  MatrixPolicy<micm::Real> other{ 4, 3, 0.0 };
  MatrixPolicy<micm::Real> other2{ 4, 3, 0.0 };

  for (micm::Index i = 0; i < 4; ++i)
  {
    for (micm::Index j = 0; j < 3; ++j)
    {
      matrix[i][j] = i * 10.3 + j * 100.5;
      other[i][j] = i * 1.7 + j * 10.2;
      other2[i][j] = i * 19.5 + j * 32.2;
    }
  }

  matrix.CopyToDevice();
  other.CopyToDevice();
  matrix.ForEach(MICM_LAMBDA(micm::Real & a, const micm::Real& b) { a += b; }, other);
  matrix.CopyToHost();
  for (micm::Index i = 0; i < 4; ++i)
  {
    for (micm::Index j = 0; j < 3; ++j)
    {
      EXPECT_NEAR(matrix[i][j], (i * 10.3 + j * 100.5) + (i * 1.7 + j * 10.2), 1.0e-5);
    }
  }

  // Reset matrix to original values
  for (micm::Index i = 0; i < 4; ++i)
  {
    for (micm::Index j = 0; j < 3; ++j)
    {
      matrix[i][j] = i * 10.3 + j * 100.5;
    }
  }

  matrix.CopyToDevice();
  other2.CopyToDevice();
  matrix.ForEach(MICM_LAMBDA(micm::Real & a, const micm::Real& b, const micm::Real& c) { a = a + b - c; }, other, other2);
  matrix.CopyToHost();
  for (micm::Index i = 0; i < 4; ++i)
  {
    for (micm::Index j = 0; j < 3; ++j)
    {
      EXPECT_NEAR(matrix[i][j], (i * 10.3 + j * 100.5) + (i * 1.7 + j * 10.2) - (i * 19.5 + j * 32.2), 1.0e-5);
    }
  }

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestSetScalar()
{
  MatrixPolicy<micm::Real> matrix{ 2, 3, 0.0 };

  matrix = 2.0;
  matrix.CopyToHost();

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 2.0);
  }

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestMax()
{
  MatrixPolicy<micm::Real> matrix{ 2, 3, 0.0 };

  matrix.Max(2.0);
  matrix.CopyToHost();

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 2.0);
  }

  matrix = 1.0;
  matrix.CopyToHost();
  matrix[1][1] = 3.0;
  matrix.CopyToDevice();
  matrix.Max(2.0);
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][0], 2.0);
  EXPECT_EQ(matrix[0][1], 2.0);
  EXPECT_EQ(matrix[0][2], 2.0);
  EXPECT_EQ(matrix[1][0], 2.0);
  EXPECT_EQ(matrix[1][1], 3.0);
  EXPECT_EQ(matrix[1][2], 2.0);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestMin()
{
  MatrixPolicy<micm::Real> matrix{ 2, 3, 0.0 };

  matrix.Min(2.0);
  matrix.CopyToHost();

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  matrix = 1.0;
  matrix.CopyToHost();
  matrix[1][1] = 3.0;
  matrix.CopyToDevice();
  matrix.Min(2.0);
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][0], 1.0);
  EXPECT_EQ(matrix[0][1], 1.0);
  EXPECT_EQ(matrix[0][2], 1.0);
  EXPECT_EQ(matrix[1][0], 1.0);
  EXPECT_EQ(matrix[1][1], 2.0);
  EXPECT_EQ(matrix[1][2], 1.0);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestPrint()
{
  MatrixPolicy<micm::Real> matrix{ 2, 3, 0.0 };

  matrix[1][1] = 3.0;

  std::stringstream ss, endline;
  ss << matrix;
  endline << std::endl;

  EXPECT_EQ(ss.str(), "0,0,0" + endline.str() + "0,3,0" + endline.str());

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestArrayFunction()
{
  using Matrix = MatrixPolicy<micm::Real>;

  Matrix matrix{ 5, 3, -1.0 };

  // Set initial values that differ by rows
  for (int i = 0; i < static_cast<int>(matrix.NumRows()); ++i)
  {
    for (int j = 0; j < static_cast<int>(matrix.NumColumns()); ++j)
    {
      matrix[i][j] = static_cast<micm::Real>(i - 2 + 10 * j);
    }
  }

  // Initial Matrix values:
  // Row 0: -2, 8, 18
  // Row 1: -1, 9, 19
  // Row 2: 0, 10, 20
  // Row 3: 1, 11, 21
  // Row 4: 2, 12, 22

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m) {
        auto tmp = m.GetRowVariable();
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, const micm::Real& c, micm::Real& t) { t = a + b + c; },
            m.GetConstColumnView(0),
            m.GetConstColumnView(1),
            m.GetConstColumnView(2),
            tmp);
        m.ForEachRow([&](micm::Real& c, const micm::Real& t) { c = 4.0 * t; }, m.GetColumnView(2), tmp);
      },
      matrix);  // pass matrix so the type and dimensions are known by the function

  matrix.CopyToDevice();
  func(matrix);  // apply the function to the matrix
  matrix.CopyToHost();

  // Check results
  EXPECT_EQ(matrix[0][2], 4.0 * (-2 + 8 + 18));  // 96
  EXPECT_EQ(matrix[1][2], 4.0 * (-1 + 9 + 19));  // 108
  EXPECT_EQ(matrix[2][2], 4.0 * (0 + 10 + 20));  // 120
  EXPECT_EQ(matrix[3][2], 4.0 * (1 + 11 + 21));  // 132
  EXPECT_EQ(matrix[4][2], 4.0 * (2 + 12 + 22));  // 144
  EXPECT_EQ(matrix[0][0], -2.0);
  EXPECT_EQ(matrix[1][0], -1.0);
  EXPECT_EQ(matrix[2][0], 0.0);
  EXPECT_EQ(matrix[3][0], 1.0);
  EXPECT_EQ(matrix[4][0], 2.0);
  EXPECT_EQ(matrix[0][1], 8.0);
  EXPECT_EQ(matrix[1][1], 9.0);
  EXPECT_EQ(matrix[2][1], 10.0);
  EXPECT_EQ(matrix[3][1], 11.0);
  EXPECT_EQ(matrix[4][1], 12.0);

  // Use the function with a different matrix with the same number of columns, but different number of rows,
  // to test that it works with different sizes
  Matrix matrix2{ 3, 3, -1.0 };
  matrix2.CopyToDevice();
  func(matrix2);
  matrix2.CopyToHost();
  EXPECT_EQ(matrix2[0][2], 4.0 * (-1 + -1 + -1));  // -12
  EXPECT_EQ(matrix2[1][2], 4.0 * (-1 + -1 + -1));  // -12
  EXPECT_EQ(matrix2[2][2], 4.0 * (-1 + -1 + -1));  // -12
  EXPECT_EQ(matrix2[0][0], -1.0);
  EXPECT_EQ(matrix2[1][0], -1.0);
  EXPECT_EQ(matrix2[2][0], -1.0);
  EXPECT_EQ(matrix2[0][1], -1.0);
  EXPECT_EQ(matrix2[1][1], -1.0);
  EXPECT_EQ(matrix2[2][1], -1.0);

  return matrix;
}

template<template<class> class MatrixPolicy>
std::tuple<MatrixPolicy<micm::Real>, MatrixPolicy<micm::Real>> TestMultiMatrixArrayFunction()
{
  using Matrix = MatrixPolicy<micm::Real>;
  Matrix matrixA{ 3, 2, 1.0 };
  Matrix matrixB{ 3, 3, 2.0 };

  // Set initial values that differ by rows
  for (micm::Index i = 0; i < matrixA.NumRows(); ++i)
  {
    for (micm::Index j = 0; j < matrixA.NumColumns(); ++j)
    {
      matrixA[i][j] = static_cast<micm::Real>(i + 10 * j);
    }
    for (micm::Index j = 0; j < matrixB.NumColumns(); ++j)
    {
      matrixB[i][j] = static_cast<micm::Real>(i * 2 + 20 * j);
    }
  }
  // Set column 2 of matrixB separately
  for (micm::Index i = 0; i < matrixB.NumRows(); ++i)
  {
    matrixB[i][2] = static_cast<micm::Real>(i * 4);
  }

  // Initial MatrixA values:
  // Row 0: 0, 10
  // Row 1: 1, 11
  // Row 2: 2, 12

  // Initial MatrixB values:
  // Row 0: 0, 20, 0
  // Row 1: 2, 22, 4
  // Row 2: 4, 24, 8

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType mA, typename Matrix::ConstViewType mB) {
        // Use an array function to set C = A + B
        // where A is from matrixA, B is from matrixB, C is in matrixA
        auto tmp = mA.GetRowVariable();
        mA.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& t) { t = a + b; },
            mA.GetConstColumnView(0),
            mB.GetConstColumnView(2),
            tmp);
        mA.ForEachRow([&](const micm::Real& t, micm::Real& c) { c = t; }, tmp, mA.GetColumnView(1));
      },
      matrixA,
      matrixB);

  matrixA.CopyToDevice();
  matrixB.CopyToDevice();
  func(matrixA, std::as_const(matrixB));
  matrixA.CopyToHost();
  matrixB.CopyToHost();

  // Check results
  EXPECT_EQ(matrixA[0][1], 0 + 0);  // 0
  EXPECT_EQ(matrixA[1][1], 1 + 4);  // 5
  EXPECT_EQ(matrixA[2][1], 2 + 8);  // 10
  EXPECT_EQ(matrixA[0][0], 0.0);
  EXPECT_EQ(matrixA[1][0], 1.0);
  EXPECT_EQ(matrixA[2][0], 2.0);
  EXPECT_EQ(matrixB[0][0], 0.0);
  EXPECT_EQ(matrixB[1][0], 2.0);
  EXPECT_EQ(matrixB[2][0], 4.0);

  return { matrixA, matrixB };
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestVectorInMatrixFunction()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;
  Matrix matrix{ 5, 3, -1.0 };

  // Set initial values that differ by rows
  for (int i = 0; i < static_cast<int>(matrix.NumRows()); ++i)
  {
    for (int j = 0; j < static_cast<int>(matrix.NumColumns()); ++j)
    {
      matrix[i][j] = static_cast<micm::Real>(i - 2 + 10 * j);
    }
  }

  // Initial Matrix values:
  // Row 0: -2, 8, 18
  // Row 1: -1, 9, 19
  // Row 2: 0, 10, 20
  // Row 3: 1, 11, 21
  // Row 4: 2, 12, 22

  // Create a vector that we will use in the function
  Vector vec(matrix.NumRows());

  // Set some initial values in the vector
  vec[0] = 100.0;
  vec[1] = 200.0;
  vec[2] = 300.0;
  vec[3] = 400.0;
  vec[4] = 500.0;

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        auto tmp = m.GetRowVariable();
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, const micm::Real& c, micm::Real& t) { t = a + b + c; },
            m.GetConstColumnView(0),
            m.GetConstColumnView(1),
            m.GetConstColumnView(2),
            tmp);
        m.ForEachRow(
            [&](micm::Real& c, const micm::Real& d, const micm::Real& t) { c = d * t; }, m.GetColumnView(2), v, tmp);
      },
      matrix,
      vec);  // pass matrix so the type and dimensions are known by the function

  matrix.CopyToDevice();
  vec.CopyToDevice();
  func(matrix, vec);  // apply the function to the matrix
  matrix.CopyToHost();

  // Check results
  EXPECT_EQ(matrix[0][2], 100.0 * (-2 + 8 + 18));  // 2400
  EXPECT_EQ(matrix[1][2], 200.0 * (-1 + 9 + 19));  // 5400
  EXPECT_EQ(matrix[2][2], 300.0 * (0 + 10 + 20));  // 9000
  EXPECT_EQ(matrix[3][2], 400.0 * (1 + 11 + 21));  // 13200
  EXPECT_EQ(matrix[4][2], 500.0 * (2 + 12 + 22));  // 17000
  EXPECT_EQ(matrix[0][0], -2.0);
  EXPECT_EQ(matrix[1][0], -1.0);
  EXPECT_EQ(matrix[2][0], 0.0);
  EXPECT_EQ(matrix[3][0], 1.0);
  EXPECT_EQ(matrix[4][0], 2.0);
  EXPECT_EQ(matrix[0][1], 8.0);
  EXPECT_EQ(matrix[1][1], 9.0);
  EXPECT_EQ(matrix[2][1], 10.0);
  EXPECT_EQ(matrix[3][1], 11.0);
  EXPECT_EQ(matrix[4][1], 12.0);

  return matrix;
}

/// @brief Test: Multiple matrices - function created with N rows, used with M rows
template<template<class> class MatrixPolicy>
std::tuple<MatrixPolicy<micm::Real>, MatrixPolicy<micm::Real>> TestMultiMatrixDifferentRowsFromCreation()
{
  using Matrix = MatrixPolicy<micm::Real>;

  // Create function with 3-row matrices
  Matrix matrixA_create{ 3, 2, 0.0 };
  Matrix matrixB_create{ 3, 3, 0.0 };

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType mA, typename Matrix::ConstViewType mB) {
        auto tmp = mA.GetRowVariable();
        mA.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& t) { t = a + b; },
            mA.GetConstColumnView(0),
            mB.GetConstColumnView(2),
            tmp);
        mA.ForEachRow([&](const micm::Real& t, micm::Real& c) { c = t * 2.0; }, tmp, mA.GetColumnView(1));
      },
      matrixA_create,
      matrixB_create);

  // Now use with 5-row matrices (different from creation)
  Matrix matrixA{ 5, 2, 0.0 };
  Matrix matrixB{ 5, 3, 0.0 };

  for (micm::Index i = 0; i < 5; ++i)
  {
    matrixA[i][0] = static_cast<micm::Real>(i + 1);
    matrixB[i][2] = static_cast<micm::Real>(i * 10);
  }

  // Should work - column counts match, row counts match each other
  matrixA.CopyToDevice();
  matrixB.CopyToDevice();
  func(matrixA, std::as_const(matrixB));
  matrixA.CopyToHost();

  EXPECT_EQ(matrixA[0][1], (1.0 + 0.0) * 2.0);   // 2
  EXPECT_EQ(matrixA[1][1], (2.0 + 10.0) * 2.0);  // 24
  EXPECT_EQ(matrixA[2][1], (3.0 + 20.0) * 2.0);  // 46
  EXPECT_EQ(matrixA[3][1], (4.0 + 30.0) * 2.0);  // 68
  EXPECT_EQ(matrixA[4][1], (5.0 + 40.0) * 2.0);  // 90

  // Also test with 2-row matrices (fewer rows than creation)
  Matrix matrixA2{ 2, 2, 0.0 };
  Matrix matrixB2{ 2, 3, 0.0 };

  matrixA2[0][0] = 10.0;
  matrixA2[1][0] = 20.0;
  matrixB2[0][2] = 5.0;
  matrixB2[1][2] = 15.0;

  matrixA2.CopyToDevice();
  matrixB2.CopyToDevice();
  func(matrixA2, std::as_const(matrixB2));
  matrixA2.CopyToHost();

  EXPECT_EQ(matrixA2[0][1], (10.0 + 5.0) * 2.0);   // 30
  EXPECT_EQ(matrixA2[1][1], (20.0 + 15.0) * 2.0);  // 70

  return { matrixA, matrixB };
}

/// @brief Test: Matrix + vector - function created with N rows, used with M rows
template<template<class> class MatrixPolicy>
std::tuple<MatrixPolicy<micm::Real>, typename MatrixPolicy<micm::Real>::template VectorType<micm::Real>>
TestMatrixVectorDifferentRowsFromCreation()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;
  // Create function with 3-row matrix and vector
  Matrix matrix_create{ 3, 3, 0.0 };
  Vector vec_create(3);

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        auto tmp = m.GetRowVariable();
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& t) { t = a + b; }, m.GetConstColumnView(0), v, tmp);
        m.ForEachRow([&](const micm::Real& t, micm::Real& c) { c = t * 3.0; }, tmp, m.GetColumnView(1));
      },
      matrix_create,
      vec_create);

  // Now use with 5-row matrix and vector (different from creation)
  Matrix matrix{ 5, 3, 0.0 };
  Vector vec(5);

  for (micm::Index i = 0; i < 5; ++i)
  {
    matrix[i][0] = static_cast<micm::Real>(i + 1);
    vec[i] = static_cast<micm::Real>(i * 10);
  }

  // Should work - columns match, row counts match each other
  matrix.CopyToDevice();
  vec.CopyToDevice();
  func(matrix, std::as_const(vec));
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][1], (1.0 + 0.0) * 3.0);   // 3
  EXPECT_EQ(matrix[1][1], (2.0 + 10.0) * 3.0);  // 36
  EXPECT_EQ(matrix[2][1], (3.0 + 20.0) * 3.0);  // 69
  EXPECT_EQ(matrix[3][1], (4.0 + 30.0) * 3.0);  // 102
  EXPECT_EQ(matrix[4][1], (5.0 + 40.0) * 3.0);  // 135

  return { matrix, vec };
}

/// @brief Test: Mismatched row counts at invocation time (should fail)
template<template<class> class MatrixPolicy>
void TestMismatchedRowsAtInvocation()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename MatrixPolicy<micm::Real>::template VectorType<micm::Real>;

  Matrix matrix_create{ 3, 2, 0.0 };
  Vector vec_create(3);

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        m.ForEachRow([&](const micm::Real& a, micm::Real& b) { b = a * 2.0; }, v, m.GetColumnView(0));
      },
      matrix_create,
      vec_create);

  // Try to invoke with matrix (5 rows) and vector (3 rows) - should fail
  Matrix matrix{ 5, 2, 0.0 };
  Vector vec(3);

  EXPECT_ANY_THROW(func(matrix, vec));

  // Try the other way - matrix (3 rows) and vector (5 rows) - should also fail
  Matrix matrix2{ 3, 2, 0.0 };
  Vector vec2(5);

  EXPECT_ANY_THROW(func(matrix2, std::as_const(vec2)));
}

/// @brief Test: Mismatched row counts between multiple matrices at invocation (should fail)
template<template<class> class MatrixPolicy>
void TestMultipleMatricesMismatchedRowsAtInvocation()
{
  using Matrix = MatrixPolicy<micm::Real>;

  Matrix matrixA_create{ 3, 2, 0.0 };
  Matrix matrixB_create{ 3, 3, 0.0 };

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType mA, typename Matrix::ConstViewType mB) {
        mA.ForEachRow(
            [&](const micm::Real& a, micm::Real& b) { b = a * 2.0; }, mB.GetConstColumnView(0), mA.GetColumnView(0));
      },
      matrixA_create,
      matrixB_create);

  // Try to invoke with matrices having different row counts - should fail
  Matrix matrixA{ 5, 2, 0.0 };
  Matrix matrixB{ 3, 3, 0.0 };  // Different row count!

  EXPECT_ANY_THROW(func(matrixA, std::as_const(matrixB)));
}

/// @brief Test: Wrong column count at invocation time (should fail)
template<template<class> class MatrixPolicy>
void TestWrongColumnCountAtInvocation()
{
  using Matrix = MatrixPolicy<micm::Real>;

  // Create function with 3-column matrix
  Matrix matrix_create{ 4, 3, 0.0 };

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m) {
        auto tmp = m.GetRowVariable();
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, const micm::Real& c, micm::Real& t) { t = a + b + c; },
            m.GetConstColumnView(0),
            m.GetConstColumnView(1),
            m.GetConstColumnView(2),
            tmp);
      },
      matrix_create);

  // Try to invoke with wrong column count - should fail
  Matrix matrix_wrong_cols{ 4, 4, 0.0 };  // 4 columns instead of 3
  EXPECT_ANY_THROW(func(matrix_wrong_cols));

  Matrix matrix_wrong_cols2{ 4, 2, 0.0 };  // 2 columns instead of 3
  EXPECT_ANY_THROW(func(matrix_wrong_cols2));

  // Should work with different row count but same column count
  Matrix matrix_ok{ 7, 3, 0.0 };  // 7 rows, 3 columns
  EXPECT_NO_THROW(func(matrix_ok));
}

template<template<class> class MatrixPolicy>
void TestMismatchedRowDimensions()
{
  using Matrix = MatrixPolicy<micm::Real>;

  Matrix matrixA{ 3, 3, 1.0 };
  Matrix matrixB{ 4, 3, 2.0 };  // Different number of rows during creation

  // Should now SUCCEED when creating with different row counts
  // (as long as column counts match, which they do here)
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType mA, typename Matrix::ConstViewType mB) {
        mA.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& c) { c = a + b; },
            mA.GetConstColumnView(0),
            mB.GetConstColumnView(0),
            mA.GetColumnView(1));
      },
      matrixA,
      matrixB);

  // Can use the function with matrices of the same row count
  Matrix matrixC{ 5, 3, 0.0 };
  Matrix matrixD{ 5, 3, 1.0 };

  EXPECT_NO_THROW(func(matrixC, std::as_const(matrixD)));

  // But should throw if invoked with mismatched row counts
  Matrix matrixE{ 3, 3, 0.0 };
  Matrix matrixF{ 4, 3, 1.0 };  // Different row count!

  EXPECT_ANY_THROW(func(matrixE, std::as_const(matrixF)));
}

template<template<class> class MatrixPolicy>
void TestMismatchedColumnDimensions()
{
  using Matrix = MatrixPolicy<micm::Real>;

  Matrix matrix{ 3, 4, 1.0 };

  // Create the function - this should succeed
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m) {
        // Try to access a column that doesn't exist
        m.ForEachRow(
            [&](const micm::Real& a, micm::Real& b) { b = a * 2.0; },
            m.GetConstColumnView(0),
            m.GetColumnView(5));  // Column 5 doesn't exist in a 4-column matrix
      },
      matrix);

  // Should fail assert when invoking the function because column 5 doesn't exist
#ifndef NDEBUG
  EXPECT_DEATH(func(matrix), "");
#endif
}

template<template<class> class MatrixPolicy>
void TestWrongMatrixDimensions()
{
  using Matrix = MatrixPolicy<micm::Real>;

  Matrix matrix1{ 3, 4, 1.0 };
  Matrix matrix2{ 3, 5, 2.0 };  // Different column count

  // Create a function that expects 4 columns
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m) {
        m.ForEachRow(
            [&](const micm::Real& a, micm::Real& b) { b = a * 2.0; },
            m.GetConstColumnView(0),
            m.GetColumnView(3));  // Column 3 exists in 4-column matrix
      },
      matrix1);

  func(matrix1);  // Should work fine
  EXPECT_NO_THROW(func(matrix1));

  // Should throw when applied to matrix with wrong column count
  EXPECT_ANY_THROW(func(matrix2));

  // But should work with different row count as long as column count matches
  Matrix matrix3{ 7, 4, 1.0 };  // 7 rows, 4 columns
  EXPECT_NO_THROW(func(matrix3));
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestMultipleTemporaries()
{
  using Matrix = MatrixPolicy<micm::Real>;

  Matrix matrix{ 4, 5, 0.0 };

  // Initialize first two columns
  for (micm::Index i = 0; i < matrix.NumRows(); ++i)
  {
    matrix[i][0] = static_cast<micm::Real>(i + 1);
    matrix[i][1] = static_cast<micm::Real>((i + 1) * 10);
  }

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m) {
        // Use TWO temporaries for intermediate calculations
        auto tmp1 = m.GetRowVariable();
        auto tmp2 = m.GetRowVariable();

        // tmp1 = col0 * col1
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& t) { t = a * b; },
            m.GetConstColumnView(0),
            m.GetConstColumnView(1),
            tmp1);

        // tmp2 = col0 + col1
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& t) { t = a + b; },
            m.GetConstColumnView(0),
            m.GetConstColumnView(1),
            tmp2);

        // col2 = tmp1 + tmp2 (product + sum)
        m.ForEachRow(
            [&](const micm::Real& t1, const micm::Real& t2, micm::Real& c) { c = t1 + t2; }, tmp1, tmp2, m.GetColumnView(2));

        // col3 = tmp1 - tmp2 (product - sum)
        m.ForEachRow(
            [&](const micm::Real& t1, const micm::Real& t2, micm::Real& c) { c = t1 - t2; }, tmp1, tmp2, m.GetColumnView(3));

        // col4 = tmp1 * tmp2
        m.ForEachRow(
            [&](const micm::Real& t1, const micm::Real& t2, micm::Real& c) { c = t1 * t2; }, tmp1, tmp2, m.GetColumnView(4));
      },
      matrix);

  matrix.CopyToDevice();
  func(matrix);
  matrix.CopyToHost();

  // Verify results
  // Row 0: col0=1, col1=10, product=10, sum=11
  EXPECT_EQ(matrix[0][2], 10.0 + 11.0);  // 21
  EXPECT_EQ(matrix[0][3], 10.0 - 11.0);  // -1
  EXPECT_EQ(matrix[0][4], 10.0 * 11.0);  // 110

  // Row 1: col0=2, col1=20, product=40, sum=22
  EXPECT_EQ(matrix[1][2], 40.0 + 22.0);  // 62
  EXPECT_EQ(matrix[1][3], 40.0 - 22.0);  // 18
  EXPECT_EQ(matrix[1][4], 40.0 * 22.0);  // 880

  // Row 3: col0=4, col1=40, product=160, sum=44
  EXPECT_EQ(matrix[3][2], 160.0 + 44.0);  // 204
  EXPECT_EQ(matrix[3][3], 160.0 - 44.0);  // 116
  EXPECT_EQ(matrix[3][4], 160.0 * 44.0);  // 7040

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestColumnViewReuse()
{
  using Matrix = MatrixPolicy<micm::Real>;

  Matrix matrix{ 3, 4, 0.0 };

  for (micm::Index i = 0; i < matrix.NumRows(); ++i)
  {
    matrix[i][0] = static_cast<micm::Real>(i + 1);
  }

  auto func = Matrix::Function(
      MICM_LAMBDA(Matrix::ViewType m) {
        // Create column views once
        auto col0 = m.GetConstColumnView(0);
        auto col1 = m.GetColumnView(1);
        auto col2 = m.GetColumnView(2);
        auto col3 = m.GetColumnView(3);

        // Reuse the same column views in multiple ForEachRow calls
        // col1 = col0 * 2
        m.ForEachRow([&](const micm::Real& a, micm::Real& b) { b = a * 2.0; }, col0, col1);

        // col2 = col0 + col1 (reusing col0 and col1)
        m.ForEachRow([&](const micm::Real& a, const micm::Real& b, micm::Real& c) { c = a + b; }, col0, col1, col2);

        // col3 = col2 * col1 (reusing col1 and col2)
        m.ForEachRow([&](const micm::Real& a, const micm::Real& b, micm::Real& c) { c = a * b; }, col2, col1, col3);
      },
      matrix);

  matrix.CopyToDevice();
  func(matrix);
  matrix.CopyToHost();

  // Row 0: col0=1, col1=2, col2=3, col3=6
  EXPECT_EQ(matrix[0][1], 2.0);
  EXPECT_EQ(matrix[0][2], 3.0);
  EXPECT_EQ(matrix[0][3], 6.0);

  // Row 1: col0=2, col1=4, col2=6, col3=24
  EXPECT_EQ(matrix[1][1], 4.0);
  EXPECT_EQ(matrix[1][2], 6.0);
  EXPECT_EQ(matrix[1][3], 24.0);

  // Row 2: col0=3, col1=6, col2=9, col3=54
  EXPECT_EQ(matrix[2][1], 6.0);
  EXPECT_EQ(matrix[2][2], 9.0);
  EXPECT_EQ(matrix[2][3], 54.0);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestFunctionReusability()
{
  using Matrix = MatrixPolicy<micm::Real>;

  Matrix matrix1{ 2, 3, 1.0 };

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m) {
        auto tmp = m.GetRowVariable();
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, const micm::Real& c, micm::Real& t) { t = a + b + c; },
            m.GetConstColumnView(0),
            m.GetConstColumnView(1),
            m.GetConstColumnView(2),
            tmp);
        m.ForEachRow([&](micm::Real& c, const micm::Real& t) { c = 2.0 * t; }, m.GetColumnView(2), tmp);
      },
      matrix1);

  // Apply to first matrix
  for (micm::Index i = 0; i < matrix1.NumRows(); ++i)
  {
    for (micm::Index j = 0; j < matrix1.NumColumns(); ++j)
    {
      matrix1[i][j] = static_cast<micm::Real>(i + j);
    }
  }

  matrix1.CopyToDevice();
  func(matrix1);
  matrix1.CopyToHost();
  EXPECT_EQ(matrix1[0][2], 2.0 * (0 + 1 + 2));  // 6
  EXPECT_EQ(matrix1[1][2], 2.0 * (1 + 2 + 3));  // 12

  // Apply to second matrix with same dimensions
  Matrix matrix2{ 2, 3, 5.0 };
  matrix2.CopyToDevice();
  func(matrix2);
  matrix2.CopyToHost();
  EXPECT_EQ(matrix2[0][2], 2.0 * (5 + 5 + 5));  // 30
  EXPECT_EQ(matrix2[1][2], 2.0 * (5 + 5 + 5));  // 30

  // Apply to third matrix with different values
  Matrix matrix3{ 2, 3, 0.0 };
  for (micm::Index i = 0; i < matrix3.NumRows(); ++i)
  {
    matrix3[i][0] = static_cast<micm::Real>(i * 10);
  }

  matrix3.CopyToDevice();
  func(matrix3);
  matrix3.CopyToHost();
  EXPECT_EQ(matrix3[0][2], 2.0 * (0 + 0 + 0));   // 0
  EXPECT_EQ(matrix3[1][2], 2.0 * (10 + 0 + 0));  // 20

  return matrix1;
}

template<template<class> class MatrixPolicy>
void TestConstMatrixFunction()
{
  using Matrix = MatrixPolicy<micm::Real>;

  Matrix matrix{ 3, 4, 0.0 };

  // Set initial values
  for (micm::Index i = 0; i < matrix.NumRows(); ++i)
  {
    for (micm::Index j = 0; j < matrix.NumColumns(); ++j)
    {
      matrix[i][j] = static_cast<micm::Real>(i * 10 + j);
    }
  }

  // Create a const reference
  const Matrix& const_matrix = matrix;

  // Create a function that only reads from the matrix
  auto read_func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ConstViewType m) {
        auto tmp = m.GetRowVariable();
        // Only use GetConstColumnView - should work with const matrices
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& t) { t = a + b; },
            m.GetConstColumnView(0),
            m.GetConstColumnView(1),
            tmp);

        // Verify we can read the values (no writes to m)
        micm::Real sum = 0.0;
        m.ForEachRow([&sum](const micm::Real& val) { sum += val; }, m.GetConstColumnView(2));
      },
      const_matrix);

  // Should work fine with const matrix
  matrix.CopyToDevice();
  EXPECT_NO_THROW(read_func(const_matrix));

  // Verify original matrix unchanged
  EXPECT_EQ(matrix[0][0], 0.0);
  EXPECT_EQ(matrix[1][2], 12.0);
}

template<template<class> class MatrixPolicy>
void TestEmptyMatrixFunction()
{
  using Matrix = MatrixPolicy<micm::Real>;

  // Test with 0 rows
  Matrix empty_rows{ 0, 3, 1.0 };

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m) {
        // This should never execute
        m.ForEachRow(
            [&](micm::Real& val) { val = 99.0; },  // Would fail if executed
            m.GetColumnView(0));
      },
      empty_rows);

  // Should not throw, just iterate 0 times
  EXPECT_NO_THROW(func(empty_rows));

  // Test with 0 columns (edge case)
  Matrix empty_cols{ 3, 0, 1.0 };

  auto func2 = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType){
          // Cannot get any column views, so just return
      },
      empty_cols);

  EXPECT_NO_THROW(func2(empty_cols));
}

// ============================================================================
// Vector Support Validation Tests
// ============================================================================

/// @brief Test: Vector with TOO FEW elements (creation succeeds, invocation fails)
template<template<class> class MatrixPolicy>
void TestVectorTooSmall()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 5, 3, 1.0 };
  Vector vec_too_small(3);  // Only 3 elements, but matrix has 5 rows

  // Should succeed at creation (row counts can differ at creation)
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ConstViewType m, typename Vector::ConstViewType v) {
        auto tmp = m.GetRowVariable();
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& t) { t = a + b; }, m.GetConstColumnView(0), v, tmp);
      },
      matrix,
      vec_too_small);

  // Should throw at invocation when vector size doesn't match matrix row count
  EXPECT_ANY_THROW(func(matrix, vec_too_small));
}

/// @brief Test: Vector with TOO MANY elements (creation succeeds, invocation fails)
template<template<class> class MatrixPolicy>
void TestVectorTooLarge()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 5, 3, 1.0 };
  Vector vec_too_large(10);  // 10 elements, but matrix has 5 rows

  // Should succeed at creation (row counts can differ at creation)
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ConstViewType m, typename Vector::ConstViewType v) {
        auto tmp = m.GetRowVariable();
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& t) { t = a + b; }, m.GetConstColumnView(0), v, tmp);
      },
      matrix,
      vec_too_large);

  // Should throw at invocation when vector size doesn't match matrix row count
  EXPECT_ANY_THROW(func(matrix, vec_too_large));
}

/// @brief Test: Empty vector with non-empty matrix (creation succeeds, invocation fails)
template<template<class> class MatrixPolicy>
void TestEmptyVectorNonEmptyMatrix()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 5, 3, 1.0 };
  Vector empty_vec;  // Empty

  // Should succeed at creation
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        m.ForEachRow([&](const micm::Real& a, micm::Real& b) { b = a; }, v, m.GetColumnView(0));
      },
      matrix,
      empty_vec);

  // Should throw at invocation when vector size doesn't match
  EXPECT_ANY_THROW(func(matrix, empty_vec));
}

/// @brief Test: Non-empty vector with empty matrix (creation succeeds, invocation fails)
template<template<class> class MatrixPolicy>
void TestNonEmptyVectorEmptyMatrix()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 0, 3, 1.0 };  // 0 rows
  Vector vec(5);

  // Should succeed at creation
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        m.ForEachRow([&](const micm::Real& a, micm::Real& b) { b = a; }, v, m.GetColumnView(0));
      },
      matrix,
      vec);

  // Should throw at invocation when vector size doesn't match matrix row count
  EXPECT_ANY_THROW(func(matrix, vec));
}

/// @brief Test: Empty vector with empty matrix (should work - no iterations)
template<template<class> class MatrixPolicy>
void TestEmptyVectorEmptyMatrix()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 0, 3, 1.0 };  // 0 rows
  Vector empty_vec;            // Empty

  // Should succeed - both are empty, ForEachRow won't iterate
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        m.ForEachRow([&](const micm::Real& a, micm::Real& b) { b = a; }, v, m.GetColumnView(0));
      },
      matrix,
      empty_vec);

  EXPECT_NO_THROW(func(matrix, empty_vec));
}

/// @brief Test: Multiple vectors with DIFFERENT sizes (creation succeeds, invocation fails)
template<template<class> class MatrixPolicy>
void TestMultipleVectorsDifferentSizes()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 5, 3, 1.0 };
  Vector vec1(5);  // Size 5
  Vector vec2(3);  // Size 3 - different!

  // Should succeed at creation (different row counts allowed at creation)
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v1, typename Vector::ConstViewType v2) {
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& c) { c = a + b; }, v1, v2, m.GetColumnView(0));
      },
      matrix,
      vec1,
      vec2);

  // Should throw at invocation because vectors have different sizes
  EXPECT_ANY_THROW(func(matrix, vec1, vec2));
}

/// @brief Test: Multiple vectors with SAME correct size (should work)
template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestMultipleVectorsSameSize()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 5, 3, 0.0 };
  Vector vec1(5);
  Vector vec2(5);

  // Initialize vectors
  for (micm::Index i = 0; i < 5; ++i)
  {
    vec1[i] = static_cast<micm::Real>(i + 1);
    vec2[i] = static_cast<micm::Real>((i + 1) * 10);
  }

  // Should succeed
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v1, typename Vector::ConstViewType v2) {
        // col0 = v1 + v2
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& c) { c = a + b; }, v1, v2, m.GetColumnView(0));
      },
      matrix,
      vec1,
      vec2);

  vec1.CopyToDevice();
  vec2.CopyToDevice();
  func(matrix, vec1, vec2);
  matrix.CopyToHost();

  // Verify results
  EXPECT_EQ(matrix[0][0], 1.0 + 10.0);  // 11
  EXPECT_EQ(matrix[1][0], 2.0 + 20.0);  // 22
  EXPECT_EQ(matrix[2][0], 3.0 + 30.0);  // 33
  EXPECT_EQ(matrix[3][0], 4.0 + 40.0);  // 44
  EXPECT_EQ(matrix[4][0], 5.0 + 50.0);  // 55

  return matrix;
}

/// @brief Test: Multiple matrices + vector - vector size must match all matrices
template<template<class> class MatrixPolicy>
std::tuple<MatrixPolicy<micm::Real>, MatrixPolicy<micm::Real>> TestMultipleMatricesOneVector()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrixA{ 4, 2, 0.0 };
  Matrix matrixB{ 4, 3, 0.0 };
  Vector vec(4);

  for (micm::Index i = 0; i < 4; ++i)
  {
    vec[i] = static_cast<micm::Real>(i * 2);
    matrixA[i][0] = static_cast<micm::Real>(i + 1);
    matrixB[i][0] = static_cast<micm::Real>(i * 10);
  }

  // Should succeed - both matrices have 4 rows, vector has 4 elements
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType mA, typename Matrix::ViewType mB, typename Vector::ConstViewType v) {
        // matrixA col1 = matrixA col0 + vector
        mA.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& c) { c = a + b; },
            mA.GetConstColumnView(0),
            v,
            mA.GetColumnView(1));

        // matrixB col1 = matrixB col0 + vector
        mB.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& c) { c = a + b; },
            mB.GetConstColumnView(0),
            v,
            mB.GetColumnView(1));
      },
      matrixA,
      matrixB,
      vec);

  matrixA.CopyToDevice();
  matrixB.CopyToDevice();
  vec.CopyToDevice();
  func(matrixA, matrixB, vec);
  matrixA.CopyToHost();
  matrixB.CopyToHost();

  // Verify results
  EXPECT_EQ(matrixA[0][1], 1.0 + 0.0);  // 1
  EXPECT_EQ(matrixA[1][1], 2.0 + 2.0);  // 4
  EXPECT_EQ(matrixA[2][1], 3.0 + 4.0);  // 7
  EXPECT_EQ(matrixA[3][1], 4.0 + 6.0);  // 10

  EXPECT_EQ(matrixB[0][1], 0.0 + 0.0);   // 0
  EXPECT_EQ(matrixB[1][1], 10.0 + 2.0);  // 12
  EXPECT_EQ(matrixB[2][1], 20.0 + 4.0);  // 24
  EXPECT_EQ(matrixB[3][1], 30.0 + 6.0);  // 36

  return { matrixA, matrixB };
}

/// @brief Test: Multiple matrices with DIFFERENT row counts + vector (creation succeeds, invocation fails)
#if defined(__GNUC__) && !defined(__clang__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
template<template<class> class MatrixPolicy>
void TestMultipleMatricesDifferentRowsVector()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrixA{ 4, 2, 0.0 };
  Matrix matrixB{ 5, 3, 0.0 };  // Different row count!
  Vector vec(4);

  // Should succeed at creation (different row counts allowed at creation)
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType mA, typename Matrix::ConstViewType mB, typename Vector::ConstViewType v) {
        mA.ForEachRow([&](const micm::Real& a, micm::Real& b) { b = a; }, v, mA.GetColumnView(0));
      },
      matrixA,
      matrixB,
      vec);

  // Should throw at invocation because matrices have different row counts
  EXPECT_ANY_THROW(func(matrixA, matrixB, vec));
}
#if defined(__GNUC__) && !defined(__clang__)
  #pragma GCC diagnostic pop
#endif

/// @brief Test: Vector size matches one matrix but not the other (creation succeeds, invocation fails)
#if defined(__GNUC__) && !defined(__clang__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
template<template<class> class MatrixPolicy>
void TestVectorSizeMatchesOneMatrixOnly()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrixA{ 5, 2, 0.0 };
  Matrix matrixB{ 5, 3, 0.0 };
  Vector vec(4);  // Wrong size for both matrices (they have 5 rows)

  // Should succeed at creation (different row counts allowed at creation)
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType mA, typename Matrix::ConstViewType mB, typename Vector::ConstViewType v) {
        mA.ForEachRow([&](const micm::Real& a, micm::Real& b) { b = a; }, v, mA.GetColumnView(0));
      },
      matrixA,
      matrixB,
      vec);

  // Should throw at invocation because vector size doesn't match matrix row counts
  EXPECT_ANY_THROW(func(matrixA, matrixB, vec));
}
#if defined(__GNUC__) && !defined(__clang__)
  #pragma GCC diagnostic pop
#endif

/// @brief Test: Const vector (read-only access)
template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestConstVector()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 3, 2, 0.0 };
  Vector vec_data = { 10.0, 20.0, 30.0 };
  const Vector& const_vec = vec_data;

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        // Read from const vector, write to matrix
        m.ForEachRow(
            [&](const micm::Real& a, micm::Real& b) { b = a * 2.0; },
            v,  // const vector
            m.GetColumnView(0));
      },
      matrix,
      const_vec);

  const_vec.CopyToDevice();
  func(matrix, const_vec);
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][0], 20.0);
  EXPECT_EQ(matrix[1][0], 40.0);
  EXPECT_EQ(matrix[2][0], 60.0);

  return matrix;
}

/// @brief Test: Non-const vector that gets modified
template<template<class> class MatrixPolicy>
std::tuple<MatrixPolicy<micm::Real>, typename MatrixPolicy<micm::Real>::template VectorType<micm::Real>> TestMutableVector()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 3, 2, 0.0 };
  Vector vec = { 5.0, 10.0, 15.0 };

  for (micm::Index i = 0; i < 3; ++i)
  {
    matrix[i][0] = static_cast<micm::Real>(i + 1);
  }

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ConstViewType m, typename Vector::ViewType v) {
        // Write to vector from matrix
        m.ForEachRow(
            [&](const micm::Real& a, micm::Real& b) { b = a * 3.0; },
            m.GetConstColumnView(0),
            v);  // non-const vector
      },
      matrix,
      vec);

  matrix.CopyToDevice();
  func(matrix, vec);
  vec.CopyToHost();

  // Vector should be modified
  EXPECT_EQ(vec[0], 3.0);
  EXPECT_EQ(vec[1], 6.0);
  EXPECT_EQ(vec[2], 9.0);

  return { matrix, vec };
}

/// @brief Test: Function reusability with different vectors
template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestFunctionReusabilityWithVectors()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 3, 2, 0.0 };
  Vector vec1 = { 1.0, 2.0, 3.0 };

  // Create function once
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        m.ForEachRow([&](const micm::Real& a, micm::Real& b) { b = a * 10.0; }, v, m.GetColumnView(0));
      },
      matrix,
      vec1);

  // Apply with first vector
  vec1.CopyToDevice();
  func(matrix, vec1);
  matrix.CopyToHost();
  EXPECT_EQ(matrix[0][0], 10.0);
  EXPECT_EQ(matrix[1][0], 20.0);
  EXPECT_EQ(matrix[2][0], 30.0);

  // Apply with second vector (same size)
  typename MatrixPolicy<micm::Real>::template VectorType<micm::Real> vec2 = { 5.0, 6.0, 7.0 };
  vec2.CopyToDevice();
  func(matrix, vec2);
  matrix.CopyToHost();
  EXPECT_EQ(matrix[0][0], 50.0);
  EXPECT_EQ(matrix[1][0], 60.0);
  EXPECT_EQ(matrix[2][0], 70.0);

  return matrix;
}

/// @brief Test: Applying function with vector of wrong size at invocation time
template<template<class> class MatrixPolicy>
void TestFunctionInvocationWithWrongSizedVector()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 3, 2, 0.0 };
  Vector vec_correct(3);
  Vector vec_wrong(5);

  // Create function with correct-sized vector
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        m.ForEachRow([&](const micm::Real& a, micm::Real& b) { b = a; }, v, m.GetColumnView(0));
      },
      matrix,
      vec_correct);

  // Should work with correct size
  EXPECT_NO_THROW(func(matrix, vec_correct));

  // Should throw when invoked with wrong-sized vector
  // This now leads to compiler warnings
  // EXPECT_ANY_THROW(func(matrix, vec_wrong));
}

/// @brief Test: Mixed - vector, column view, and row variable together
template<template<class> class MatrixPolicy>
MatrixPolicy<micm::Real> TestMixedVectorColumnViewRowVariable()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 4, 3, 0.0 };
  Vector vec(4);

  for (micm::Index i = 0; i < 4; ++i)
  {
    matrix[i][0] = static_cast<micm::Real>(i + 1);
    matrix[i][1] = static_cast<micm::Real>((i + 1) * 10);
    vec[i] = static_cast<micm::Real>((i + 1) * 100);
  }

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        auto tmp = m.GetRowVariable();

        // tmp = col0 + col1
        m.ForEachRow(
            [&](const micm::Real& a, const micm::Real& b, micm::Real& t) { t = a + b; },
            m.GetConstColumnView(0),
            m.GetConstColumnView(1),
            tmp);

        // col2 = tmp + vector
        m.ForEachRow(
            [&](const micm::Real& t, const micm::Real& v_elem, micm::Real& result) { result = t + v_elem; },
            tmp,
            v,
            m.GetColumnView(2));
      },
      matrix,
      vec);

  matrix.CopyToDevice();
  vec.CopyToDevice();
  func(matrix, vec);
  matrix.CopyToHost();
  vec.CopyToHost();

  // Row 0: col0=1, col1=10, vec=100, result=(1+10)+100=111
  EXPECT_EQ(matrix[0][2], 111.0);
  // Row 1: col0=2, col1=20, vec=200, result=(2+20)+200=222
  EXPECT_EQ(matrix[1][2], 222.0);
  // Row 2: col0=3, col1=30, vec=300, result=(3+30)+300=333
  EXPECT_EQ(matrix[2][2], 333.0);
  // Row 3: col0=4, col1=40, vec=400, result=(4+40)+400=444
  EXPECT_EQ(matrix[3][2], 444.0);

  return matrix;
}

/// @brief Test: Different integer types (typename MatrixPolicy<micm::Real>::template VectorType<int>)
template<template<class> class MatrixPolicy>
MatrixPolicy<int> TestIntegerVector()
{
  using Matrix = MatrixPolicy<int>;
  using Vector = typename Matrix::template VectorType<int>;

  Matrix matrix{ 3, 2, 0 };
  Vector vec = { 10, 20, 30 };

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        m.ForEachRow([&](const int& a, int& b) { b = a * 2; }, v, m.GetColumnView(0));
      },
      matrix,
      vec);

  vec.CopyToDevice();
  func(matrix, vec);
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][0], 20);
  EXPECT_EQ(matrix[1][0], 40);
  EXPECT_EQ(matrix[2][0], 60);

  return matrix;
}

template<template<class> class MatrixPolicy>
void TestFunctionWithConstSignature()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 3, 2, 0.0 };
  Vector vec(3);

  // Create function
  auto func_auto = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) {
        m.ForEachRow([&](const micm::Real& a, micm::Real& b) { b = a * 2.0; }, v, m.GetColumnView(0));
      },
      matrix,
      vec);

  // Try to wrap in std::function with const signature
  std::function<void(Matrix&, const Vector&)> func_std = func_auto;

  func_std(matrix, vec);
}

template<template<class> class MatrixPolicy>
void TestFill()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 3, 2, 0.0 };

  // Fill a matrix column with a scalar value.
  {
    auto func = Matrix::Function(MICM_LAMBDA(typename Matrix::ViewType m) { m.Fill(m.GetColumnView(1), 3.2); }, matrix);

    func(matrix);
    matrix.CopyToHost();

    EXPECT_EQ(matrix[0][0], 0.0);
    EXPECT_EQ(matrix[0][1], static_cast<micm::Real>(3.2));
    EXPECT_EQ(matrix[1][0], 0.0);
    EXPECT_EQ(matrix[1][1], static_cast<micm::Real>(3.2));
    EXPECT_EQ(matrix[2][0], 0.0);
    EXPECT_EQ(matrix[2][1], static_cast<micm::Real>(3.2));
  }

  // Fill a caller-owned typename MatrixPolicy<micm::Real>::template VectorType with a scalar value.
  {
    Vector vec(3);
    auto func = Matrix::Function(
        MICM_LAMBDA(typename Matrix::ConstViewType m, typename Vector::ViewType v) { m.Fill(v, 3.2); }, matrix, vec);

    func(matrix, vec);
    vec.CopyToHost();

    EXPECT_EQ(vec[0], static_cast<micm::Real>(3.2));
    EXPECT_EQ(vec[1], static_cast<micm::Real>(3.2));
    EXPECT_EQ(vec[2], static_cast<micm::Real>(3.2));
  }

  // Fill a caller-owned row-variable temp with a scalar value.
  {
    auto func = Matrix::Function(
        MICM_LAMBDA(typename Matrix::ViewType m) {
          auto tmp = m.GetRowVariable();
          m.Fill(tmp, 9.9);
          // Broadcast the temp into column 0 so we can observe it from outside.
          m.ForEachRow([](micm::Real& c, const micm::Real& t) { c = t; }, m.GetColumnView(0), tmp);
        },
        matrix);

    func(matrix);
    matrix.CopyToHost();

    EXPECT_EQ(matrix[0][0], static_cast<micm::Real>(9.9));
    EXPECT_EQ(matrix[1][0], static_cast<micm::Real>(9.9));
    EXPECT_EQ(matrix[2][0], static_cast<micm::Real>(9.9));
  }
}

template<template<class> class MatrixPolicy>
void TestCopy()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Real>;

  Matrix matrix{ 3, 2, 0.0 };
  Vector vec{ 3.2, 4.2, 1.3 };

  // Copy a typename MatrixPolicy<micm::Real>::template VectorType into a matrix column.
  {
    vec.CopyToDevice();
    matrix.CopyToDevice();
    auto func = Matrix::Function(
        MICM_LAMBDA(typename Matrix::ViewType m, typename Vector::ConstViewType v) { m.Copy(m.GetColumnView(1), v); },
        matrix,
        vec);

    func(matrix, vec);
    matrix.CopyToHost();

    EXPECT_EQ(matrix[0][0], 0.0);
    EXPECT_EQ(matrix[0][1], static_cast<micm::Real>(3.2));
    EXPECT_EQ(matrix[1][0], 0.0);
    EXPECT_EQ(matrix[1][1], static_cast<micm::Real>(4.2));
    EXPECT_EQ(matrix[2][0], 0.0);
    EXPECT_EQ(matrix[2][1], static_cast<micm::Real>(1.3));
  }

  // Copy a const matrix column into a typename MatrixPolicy<micm::Real>::template VectorType.
  {
    Vector vec2(3);
    vec2.CopyToDevice();
    auto func = Matrix::Function(
        MICM_LAMBDA(typename Matrix::ConstViewType m, typename Vector::ViewType v) { m.Copy(v, m.GetConstColumnView(1)); },
        matrix,
        vec2);

    func(matrix, vec2);
    vec2.CopyToHost();

    EXPECT_EQ(vec2[0], static_cast<micm::Real>(3.2));
    EXPECT_EQ(vec2[1], static_cast<micm::Real>(4.2));
    EXPECT_EQ(vec2[2], static_cast<micm::Real>(1.3));
  }

  // Copy one matrix column into another (mutable-to-mutable).
  {
    auto func = Matrix::Function(
        MICM_LAMBDA(typename Matrix::ViewType m) { m.Copy(m.GetColumnView(0), m.GetColumnView(1)); }, matrix);

    func(matrix);
    matrix.CopyToHost();

    EXPECT_EQ(matrix[0][0], static_cast<micm::Real>(3.2));
    EXPECT_EQ(matrix[0][1], static_cast<micm::Real>(3.2));
    EXPECT_EQ(matrix[1][0], static_cast<micm::Real>(4.2));
    EXPECT_EQ(matrix[1][1], static_cast<micm::Real>(4.2));
    EXPECT_EQ(matrix[2][0], static_cast<micm::Real>(1.3));
    EXPECT_EQ(matrix[2][1], static_cast<micm::Real>(1.3));
  }

  // Copy one matrix column into another (const-to-mutable via GetConstColumnView).
  {
    // Reset column 0 so we can observe the copy.
    for (std::size_t i = 0; i < matrix.NumRows(); ++i)
    {
      matrix[i][0] = 0.0;
    }

    auto func = Matrix::Function(
        MICM_LAMBDA(typename Matrix::ViewType m) { m.Copy(m.GetColumnView(0), m.GetConstColumnView(1)); }, matrix);

    func(matrix);
    matrix.CopyToHost();

    EXPECT_EQ(matrix[0][0], static_cast<micm::Real>(3.2));
    EXPECT_EQ(matrix[1][0], static_cast<micm::Real>(4.2));
    EXPECT_EQ(matrix[2][0], static_cast<micm::Real>(1.3));
  }

  // Round-trip: matrix column -> row-variable temp -> matrix column.
  {
    auto func = Matrix::Function(
        MICM_LAMBDA(typename Matrix::ViewType m) {
          auto tmp = m.GetRowVariable();
          m.Copy(tmp, m.GetConstColumnView(1));
          // Zero column 1 first so the copy-back is observable.
          m.Fill(m.GetColumnView(1), 0.0);
          m.ForEachRow([](micm::Real& c, const micm::Real& t) { c = t; }, m.GetColumnView(1), tmp);
        },
        matrix);

    func(matrix);
    matrix.CopyToHost();

    EXPECT_EQ(matrix[0][1], static_cast<micm::Real>(3.2));
    EXPECT_EQ(matrix[1][1], static_cast<micm::Real>(4.2));
    EXPECT_EQ(matrix[2][1], static_cast<micm::Real>(1.3));
  }
}

/// @brief Reduce (Sum): sum of x² across every element of the matrix, driven through
///        Function()'s GroupView. Verifies that non-strict Reduce (which touches
///        any trailing padding cells for L>1) still produces the correct result
///        because padding cells were initialized to the reducer's identity (0).
template<template<class> class MatrixPolicy>
void TestReduceSum()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Scalar = typename Matrix::template ScalarType<micm::Real>;
  using Sum = typename Matrix::template SumType<micm::Real>;

  Matrix matrix{ 3, 2, 0.0 };
  matrix[0][0] = 1.0;
  matrix[0][1] = 2.0;
  matrix[1][0] = 3.0;
  matrix[1][1] = 4.0;
  matrix[2][0] = 5.0;
  matrix[2][1] = 6.0;
  matrix.CopyToDevice();

  Scalar total = 0.0;
  total.CopyToDevice();
  Sum total_sum{ total };  // must construct outside of lambda (on host)
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ConstViewType view) {
        view.Reduce(total_sum, [](const micm::Real& a, micm::Real& acc) { acc += a * a; }, view.GetConstColumnView(0));
        view.Reduce(total_sum, [](const micm::Real& a, micm::Real& acc) { acc += a * a; }, view.GetConstColumnView(1));
      },
      matrix);
  func(matrix);
  total.CopyToHost();

  // 1 + 4 + 9 + 16 + 25 + 36 = 91
  EXPECT_EQ(total, 91.0);
}

/// @brief Reduce (Max): max element across every column of the matrix.
template<template<class> class MatrixPolicy>
void TestReduceMax()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Scalar = typename Matrix::template ScalarType<micm::Real>;
  using Max = typename Matrix::template MaxType<micm::Real>;

  Matrix matrix{ 3, 2, 0.0 };
  matrix[0][0] = 1.0;
  matrix[0][1] = 8.0;
  matrix[1][0] = 5.0;
  matrix[1][1] = 2.0;
  matrix[2][0] = 3.0;
  matrix[2][1] = 4.0;
  matrix.CopyToDevice();

  Scalar max_val = std::numeric_limits<micm::Real>::lowest();
  max_val.CopyToDevice();
  Max max_val_max{ max_val };  // must construct outside of lambda
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ConstViewType view) {
        view.Reduce(
            max_val_max,
            [](const micm::Real& a, micm::Real& acc)
            {
              if (a > acc)
              {
                acc = a;
              }
            },
            view.GetConstColumnView(0));
        view.Reduce(
            max_val_max,
            [](const micm::Real& a, micm::Real& acc)
            {
              if (a > acc)
              {
                acc = a;
              }
            },
            view.GetConstColumnView(1));
      },
      matrix);
  func(matrix);
  max_val.CopyToHost();

  EXPECT_EQ(max_val, 8.0);
}

/// @brief Reduce (LOr): true iff any element is greater than a threshold.
template<template<class> class MatrixPolicy>
void TestReduceLOr()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Scalar = typename Matrix::template ScalarType<micm::Bool>;
  using LOr = typename Matrix::LOrType;

  // Case 1: no element exceeds the threshold -> LOr result stays false.
  {
    Matrix matrix{ 3, 2, 1.0 };
    matrix.CopyToDevice();

    Scalar any_large = false;
    any_large.CopyToDevice();
    LOr any_large_lor{ any_large };  // must construct outside of lambda
    auto func = Matrix::Function(
        MICM_LAMBDA(typename Matrix::ConstViewType view) {
          view.Reduce(
              any_large_lor,
              [](const micm::Real& a, micm::Bool& acc) { acc = acc || (a > 10.0); },
              view.GetConstColumnView(0));
          view.Reduce(
              any_large_lor,
              [](const micm::Real& a, micm::Bool& acc) { acc = acc || (a > 10.0); },
              view.GetConstColumnView(1));
        },
        matrix);
    func(matrix);
    any_large.CopyToHost();

    EXPECT_FALSE(any_large);
  }

  // Case 2: one element exceeds the threshold -> LOr result becomes true.
  {
    Matrix matrix{ 3, 2, 1.0 };
    matrix.CopyToHost();
    matrix[2][1] = 42.0;
    matrix.CopyToDevice();

    Scalar any_large = false;
    any_large.CopyToDevice();
    LOr any_large_lor{ any_large };  // must construct outside of lambda
    auto func = Matrix::Function(
        MICM_LAMBDA(typename Matrix::ConstViewType view) {
          view.Reduce(
              any_large_lor,
              [](const micm::Real& a, micm::Bool& acc) { acc = acc || (a > 10.0); },
              view.GetConstColumnView(0));
          view.Reduce(
              any_large_lor,
              [](const micm::Real& a, micm::Bool& acc) { acc = acc || (a > 10.0); },
              view.GetConstColumnView(1));
        },
        matrix);
    func(matrix);
    any_large.CopyToHost();

    EXPECT_TRUE(any_large);
  }
}

/// @brief Reduce (LAnd): true if every element satisfies a predicate.
///        The predicate `std::isfinite(x)` also holds for the (zero-initialized)
///        padding cells that a non-strict Reduce may touch on L>1 policies, so
///        this can use Reduce (not ReduceStrict).
template<template<class> class MatrixPolicy>
void TestReduceLAnd()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Scalar = typename Matrix::template ScalarType<micm::Bool>;
  using LAnd = typename Matrix::LAndType;

  // Case 1: all elements finite -> LAnd result stays true.
  {
    Matrix matrix{ 3, 2, 1.0 };
    matrix.CopyToDevice();

    Scalar all_finite = true;
    all_finite.CopyToDevice();
    LAnd all_finite_land{ all_finite };  // must construct outside of lambda (on host)
    auto func = Matrix::Function(
        MICM_LAMBDA(typename Matrix::ConstViewType view) {
          view.Reduce(
              all_finite_land,
              [](const micm::Real& a, micm::Bool& acc) { acc = acc && std::isfinite(a); },
              view.GetConstColumnView(0));
          view.Reduce(
              all_finite_land,
              [](const micm::Real& a, micm::Bool& acc) { acc = acc && std::isfinite(a); },
              view.GetConstColumnView(1));
        },
        matrix);
    func(matrix);
    all_finite.CopyToHost();

    EXPECT_TRUE(all_finite);
  }

  // Case 2: one element non-finite -> LAnd result becomes false.
  {
    Matrix matrix{ 3, 2, 1.0 };
    matrix.CopyToHost();
    matrix[1][0] = std::numeric_limits<micm::Real>::quiet_NaN();
    matrix.CopyToDevice();

    Scalar all_finite = true;
    all_finite.CopyToDevice();
    LAnd all_finite_land{ all_finite };  // must construct outside of lambda (on host)
    auto func = Matrix::Function(
        MICM_LAMBDA(typename Matrix::ConstViewType view) {
          view.Reduce(
              all_finite_land,
              [](const micm::Real& a, micm::Bool& acc) { acc = acc && std::isfinite(a); },
              view.GetConstColumnView(0));
          view.Reduce(
              all_finite_land,
              [](const micm::Real& a, micm::Bool& acc) { acc = acc && std::isfinite(a); },
              view.GetConstColumnView(1));
        },
        matrix);
    func(matrix);
    all_finite.CopyToHost();

    EXPECT_FALSE(all_finite);
  }
}

/// @brief ReduceStrict: verify the strict variant only visits real (non-padding)
///        rows even when L doesn't evenly divide the row count.
template<template<class> class MatrixPolicy>
void TestReduceStrict()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Scalar = typename Matrix::template ScalarType<micm::Real>;
  using Sum = typename Matrix::template SumType<micm::Real>;

  // 3 rows, but for L > 1 policies the tail group has < L real rows plus padding.
  // ReduceStrict with a counter-like lambda should visit exactly NumRows() real
  // rows -- never any padding rows.
  Matrix matrix{ 3, 1, 0.0 };
  matrix.CopyToDevice();

  Scalar count = 0.0;
  count.CopyToDevice();
  Sum count_sum{ count };  // must construct outside of lambda (on host)
  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ConstViewType view) {
        view.ReduceStrict(count_sum, [](const micm::Real&, micm::Real& acc) { acc += 1.0; }, view.GetConstColumnView(0));
      },
      matrix);
  func(matrix);
  count.CopyToHost();

  EXPECT_EQ(count, 3.0);
}

/// @brief Vector capture: verify that a vector view can be captured and
///        used in Function() lambdas
template<template<class> class MatrixPolicy>
void TestVectorCapture()
{
  using Matrix = MatrixPolicy<micm::Real>;
  using Vector = typename Matrix::template VectorType<micm::Index>;

  Matrix matrix{ 3, 3, 0.0 };
  matrix.CopyToDevice();

  Vector vec1{ 2, 0, 1 };
  Vector vec2{ 10, 20, 30 };
  vec1.CopyToDevice();
  vec2.CopyToDevice();

  auto vec1_view = vec1.GetView();
  auto vec2_view = std::as_const(vec2).GetView();

  auto func = Matrix::Function(
      MICM_LAMBDA(typename Matrix::ViewType m) {
        for (auto v1 : vec1_view)
        {
          m.ForEachRow([=](micm::Real& a) { a += vec2_view[v1]; }, m.GetColumnView(0));
        }
        m.ForEachRow([=](micm::Real& a) { a += vec2_view[vec1_view[1]]; }, m.GetColumnView(1));
      },
      matrix);
  func(matrix);
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][0], 60.0);
  EXPECT_EQ(matrix[0][1], 10.0);
  EXPECT_EQ(matrix[0][2], 0.0);
  EXPECT_EQ(matrix[1][0], 60.0);
  EXPECT_EQ(matrix[1][1], 10.0);
  EXPECT_EQ(matrix[1][2], 0.0);
  EXPECT_EQ(matrix[2][0], 60.0);
  EXPECT_EQ(matrix[2][1], 10.0);
  EXPECT_EQ(matrix[2][2], 0.0);
}
