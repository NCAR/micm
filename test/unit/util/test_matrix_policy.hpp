// Tests of common matrix functions
#include <gtest/gtest.h>

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testSmallMatrix()
{
  MatrixPolicy<double> matrix(3, 5);

  matrix[1][3] = 64.7;
  matrix[0][0] = 41.2;
  matrix[2][4] = 102.3;

  EXPECT_EQ(matrix[1][3], 64.7);
  EXPECT_EQ(matrix[0][0], 41.2);
  EXPECT_EQ(matrix[2][4], 102.3);

  std::vector<double> &data = matrix.AsVector();

  EXPECT_GE(data.size(), 3);

  return matrix;
}

template<template<class> class MatrixPolicy>
const MatrixPolicy<double> testSmallConstMatrix()
{
  MatrixPolicy<double> matrix(3, 5);

  matrix[1][3] = 64.7;
  matrix[0][0] = 41.2;
  matrix[2][4] = 102.3;

  const MatrixPolicy<double> const_matrix = matrix;

  EXPECT_EQ(const_matrix[1][3], 64.7);
  EXPECT_EQ(const_matrix[0][0], 41.2);
  EXPECT_EQ(const_matrix[2][4], 102.3);

  const std::vector<double> &data = const_matrix.AsVector();

  EXPECT_GE(data.size(), 3);

  return const_matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testInializeMatrix()
{
  MatrixPolicy<double> matrix{ 2, 3, 12.4 };

  EXPECT_EQ(matrix[0][0], 12.4);
  EXPECT_EQ(matrix[1][0], 12.4);
  EXPECT_EQ(matrix[1][2], 12.4);

  return matrix;
}

template<template<class> class MatrixPolicy>
const MatrixPolicy<double> testInializeConstMatrix()
{
  const MatrixPolicy<double> matrix{ 2, 3, 12.4 };

  EXPECT_EQ(matrix[0][0], 12.4);
  EXPECT_EQ(matrix[1][0], 12.4);
  EXPECT_EQ(matrix[1][2], 12.4);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<int> testLoopOverMatrix()
{
  MatrixPolicy<int> matrix(3, 4, 0);
  for (std::size_t i{}; i < matrix.NumRows(); ++i)
  {
    for (std::size_t j{}; j < matrix.NumColumns(); ++j)
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
const MatrixPolicy<int> testLoopOverConstMatrix()
{
  MatrixPolicy<int> matrix(3, 4, 0);
  for (std::size_t i{}; i < matrix.NumRows(); ++i)
  {
    for (std::size_t j{}; j < matrix.NumColumns(); ++j)
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
MatrixPolicy<int> testStrides()
{
  MatrixPolicy<int> matrix(3, 4, 0);

  for (std::size_t i = 0; i < matrix.NumRows(); ++i)
    for (std::size_t j = 0; j < matrix.NumColumns(); ++j)
      matrix.AsVector()[i * matrix.RowStride() + j * matrix.ColumnStride()] = i * 100 + j;

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
MatrixPolicy<double> testConversionToVector()
{
  MatrixPolicy<double> matrix{ 2, 3, 0.0 };

  matrix[1][0] = 13.2;
  matrix[1][1] = 31.2;
  matrix[1][2] = 314.2;

  std::vector<double> slice = matrix[1];

  EXPECT_EQ(slice[0], 13.2);
  EXPECT_EQ(slice[1], 31.2);
  EXPECT_EQ(slice[2], 314.2);

  return matrix;
}

template<template<class> class MatrixPolicy>
const MatrixPolicy<double> testConstConversionToVector()
{
  MatrixPolicy<double> matrix{ 2, 3, 0.0 };

  matrix[1][0] = 13.2;
  matrix[1][1] = 31.2;
  matrix[1][2] = 314.2;

  const MatrixPolicy<double> const_matrix = matrix;
  std::vector<double> slice = const_matrix[1];

  EXPECT_EQ(slice[0], 13.2);
  EXPECT_EQ(slice[1], 31.2);
  EXPECT_EQ(slice[2], 314.2);

  return const_matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testConversionFromVector()
{
  MatrixPolicy<double> zero_matrix = std::vector<std::vector<double>>{};

  EXPECT_EQ(zero_matrix.NumRows(), 0);

  std::vector<std::vector<double>> vec = { { 412.3, 32.4, 41.3 }, { 5.33, -0.3, 31.2 } };

  MatrixPolicy<double> matrix = vec;

  EXPECT_EQ(matrix.NumRows(), 2);
  EXPECT_EQ(matrix.NumColumns(), 3);
  EXPECT_EQ(matrix[0].Size(), 3);
  EXPECT_EQ(matrix[0][0], 412.3);
  EXPECT_EQ(matrix[0][1], 32.4);
  EXPECT_EQ(matrix[0][2], 41.3);
  EXPECT_EQ(matrix[1].Size(), 3);
  EXPECT_EQ(matrix[1][0], 5.33);
  EXPECT_EQ(matrix[1][1], -0.3);
  EXPECT_EQ(matrix[1][2], 31.2);

  std::vector<std::vector<int>> bad_vector = { { 3 }, { 4, 5 }, { 5 } };

  MatrixPolicy<int> bad_matrix;
  EXPECT_ANY_THROW(bad_matrix = bad_vector);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testAssignmentFromVector()
{
  MatrixPolicy<double> matrix{ 4, 3, 0.0 };
  std::vector<double> other = { 12.3, 15.1, 24.3 };
  std::vector<double> big_other = { 14.3, 52.3, 65.7, 16.34 };
  std::vector<double> small_other = { 13.2, 52.8 };

  matrix[2] = other;

  EXPECT_EQ(matrix[0][0], 0.0);
  EXPECT_EQ(matrix[2][0], 12.3);
  EXPECT_EQ(matrix[2][1], 15.1);
  EXPECT_EQ(matrix[2][2], 24.3);
  EXPECT_EQ(matrix[3][0], 0.0);

  matrix[2] = big_other;

  EXPECT_EQ(matrix[0][0], 0.0);
  EXPECT_EQ(matrix[2][0], 14.3);
  EXPECT_EQ(matrix[2][1], 52.3);
  EXPECT_EQ(matrix[2][2], 65.7);
  EXPECT_EQ(matrix[3][0], 0.0);

  EXPECT_ANY_THROW(matrix[2] = small_other);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testAxpy()
{
  std::size_t num_rows = 4;
  std::size_t num_columns = 3;
  MatrixPolicy<double> matrix{ num_rows, num_columns, 100.0 };
  MatrixPolicy<double> other{ num_rows, num_columns, 200.0 };
  double alpha = 1.39;
  double sum = 0.0;
  double result = 0.0;

  for (int i = 0; i < num_rows; ++i)
    for (int j = 0; j < num_columns; ++j)
    {
      auto y = i * 10.3 + j * 100.5;
      auto x = i * 1.7 + j * 10.2;
      matrix[i][j] = y;
      other[i][j] = x;
      sum += y + alpha * x;
    }

  matrix.Axpy(alpha, other);

  for (int i = 0; i < num_rows; ++i)
    for (int j = 0; j < num_columns; ++j)
      result += matrix[i][j];
  EXPECT_NEAR(sum, result, 1.0e-5);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testForEach()
{
  MatrixPolicy<double> matrix{ 4, 3, 100.0 };
  MatrixPolicy<double> other{ 4, 3, 200.0 };
  MatrixPolicy<double> other2{ 4, 3, 300.0 };
  double sum = 0.0;
  double sum2 = 0.0;
  double result = 0.0;

  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 3; ++j)
    {
      matrix[i][j] = i * 10.3 + j * 100.5;
      other[i][j] = i * 1.7 + j * 10.2;
      sum += i * 10.3 + j * 100.5 + i * 1.7 + j * 10.2;
      other2[i][j] = i * 19.5 + j * 32.2;
      sum2 += i * 10.3 + j * 100.5 + i * 1.7 + j * 10.2 - i * 19.5 - j * 32.2;
    }

  matrix.ForEach([&](double &a, const double &b) { result += a + b; }, other);
  EXPECT_NEAR(sum, result, 1.0e-5);
  result = 0.0;
  matrix.ForEach([&](double &a, const double &b, const double &c) { result += a + b - c; }, other, other2);
  EXPECT_NEAR(sum2, result, 1.0e-5);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testSetScalar()
{
  MatrixPolicy<double> matrix{ 2, 3, 0.0 };

  matrix = 2.0;

  for (auto &elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 2.0);
  }

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testMax()
{
  MatrixPolicy<double> matrix{ 2, 3, 0.0 };

  matrix.Max(2.0);

  for (auto &elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 2.0);
  }

  matrix = 1.0;
  matrix[1][1] = 3.0;
  matrix.Max(2.0);

  EXPECT_EQ(matrix[0][0], 2.0);
  EXPECT_EQ(matrix[0][1], 2.0);
  EXPECT_EQ(matrix[0][2], 2.0);
  EXPECT_EQ(matrix[1][0], 2.0);
  EXPECT_EQ(matrix[1][1], 3.0);
  EXPECT_EQ(matrix[1][2], 2.0);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testMin()
{
  MatrixPolicy<double> matrix{ 2, 3, 0.0 };

  matrix.Min(2.0);

  for (auto &elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  matrix = 1.0;
  matrix[1][1] = 3.0;
  matrix.Min(2.0);

  EXPECT_EQ(matrix[0][0], 1.0);
  EXPECT_EQ(matrix[0][1], 1.0);
  EXPECT_EQ(matrix[0][2], 1.0);
  EXPECT_EQ(matrix[1][0], 1.0);
  EXPECT_EQ(matrix[1][1], 2.0);
  EXPECT_EQ(matrix[1][2], 1.0);

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testPrint()
{
  MatrixPolicy<double> matrix{ 2, 3, 0.0 };

  matrix[1][1] = 3.0;

  std::stringstream ss, endline;
  ss << matrix;
  endline << std::endl;

  EXPECT_EQ(ss.str(), "0,0,0" + endline.str() + "0,3,0" + endline.str());

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testArrayFunction()
{
  MatrixPolicy<double> matrix{ 5, 3, -1.0 };

  // Set initial values that differ by rows
  for (int i = 0; i < static_cast<int>(matrix.NumRows()); ++i)
    for (int j = 0; j < static_cast<int>(matrix.NumColumns()); ++j)
      matrix[i][j] = static_cast<double>(i - 2 + 10 * j);

  // Initial Matrix values:
  // Row 0: -2, 8, 18
  // Row 1: -1, 9, 19
  // Row 2: 0, 10, 20
  // Row 3: 1, 11, 21
  // Row 4: 2, 12, 22

  auto func = MatrixPolicy<double>::Function(
    [](auto&& m)
    {
      auto tmp = m.GetRowVariable();
      m.ForEachRow([&](const double& a, const double& b, const double& c, double& t)
        { t = a + b + c; },
        m.GetConstColumnView(0),
        m.GetConstColumnView(1),
        m.GetConstColumnView(2),
        tmp);
      m.ForEachRow([&](double& c, const double& t)
        { c = 4.0 * t; },
        m.GetColumnView(2),
        tmp);
    }, matrix); // pass matrix so the type and dimensions are known by the function

  func(matrix); // apply the function to the matrix

  // Check results
  EXPECT_EQ(matrix[0][2], 4.0 * (-2 + 8 + 18));   // 96
  EXPECT_EQ(matrix[1][2], 4.0 * (-1 + 9 + 19));   // 108
  EXPECT_EQ(matrix[2][2], 4.0 * (0 + 10 + 20));   // 120
  EXPECT_EQ(matrix[3][2], 4.0 * (1 + 11 + 21));   // 132
  EXPECT_EQ(matrix[4][2], 4.0 * (2 + 12 + 22));   // 144
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

  // Use the function with a different matrix with the same dimensions
  MatrixPolicy<double> matrix2{ 5, 3, -1.0 };
  func(matrix2);
  EXPECT_EQ(matrix2[0][2], 4.0 * (-1 + -1 + -1));   // -12
  EXPECT_EQ(matrix2[1][2], 4.0 * (-1 + -1 + -1));   // -12
  EXPECT_EQ(matrix2[2][2], 4.0 * (-1 + -1 + -1));   // -12
  EXPECT_EQ(matrix2[3][2], 4.0 * (-1 + -1 + -1));   // -12
  EXPECT_EQ(matrix2[4][2], 4.0 * (-1 + -1 + -1));   // -12
  EXPECT_EQ(matrix2[0][0], -1.0);
  EXPECT_EQ(matrix2[1][0], -1.0);
  EXPECT_EQ(matrix2[2][0], -1.0);
  EXPECT_EQ(matrix2[3][0], -1.0);
  EXPECT_EQ(matrix2[4][0], -1.0);
  EXPECT_EQ(matrix2[0][1], -1.0);
  EXPECT_EQ(matrix2[1][1], -1.0);
  EXPECT_EQ(matrix2[2][1], -1.0);
  EXPECT_EQ(matrix2[3][1], -1.0);
  EXPECT_EQ(matrix2[4][1], -1.0);

  return matrix;
}

template<template<class> class MatrixPolicy>
std::tuple<MatrixPolicy<double>, MatrixPolicy<double>> testMultiMatrixArrayFunction()
{
  MatrixPolicy<double> matrixA{ 3, 2, 1.0 };
  MatrixPolicy<double> matrixB{ 3, 3, 2.0 };

  // Set initial values that differ by rows
  for (int i = 0; i < static_cast<int>(matrixA.NumRows()); ++i)
  {
    for (int j = 0; j < static_cast<int>(matrixA.NumColumns()); ++j)
    {
      matrixA[i][j] = static_cast<double>(i + 10 * j);
    }
    for (int j = 0; j < static_cast<int>(matrixB.NumColumns()); ++j)
    {
      matrixB[i][j] = static_cast<double>(i * 2 + 20 * j);
    }
  }
  // Set column 2 of matrixB separately
  for (int i = 0; i < static_cast<int>(matrixB.NumRows()); ++i)
    matrixB[i][2] = static_cast<double>(i * 4);

  // Initial MatrixA values:
  // Row 0: 0, 10
  // Row 1: 1, 11
  // Row 2: 2, 12

  // Initial MatrixB values:
  // Row 0: 0, 20, 0
  // Row 1: 2, 22, 4
  // Row 2: 4, 24, 8

  auto func = MatrixPolicy<double>::Function(
    [](auto&& mA, auto&& mB)
    {
      // Use an array function to set C = A + B
      // where A is from matrixA, B is from matrixB, C is in matrixA
      auto tmp = mA.GetRowVariable();
      mA.ForEachRow([&](const double& a, const double& b, double& t)
        { t = a + b; },
        mA.GetConstColumnView(0),
        mB.GetConstColumnView(2),
        tmp);
      mA.ForEachRow([&](const double& t, double& c)
        { c = t; },
        tmp,
        mA.GetColumnView(1));
    }, matrixA, matrixB);

  func(matrixA, matrixB);

  // Check results
  EXPECT_EQ(matrixA[0][1], 0 + 0);   // 0
  EXPECT_EQ(matrixA[1][1], 1 + 4);   // 5
  EXPECT_EQ(matrixA[2][1], 2 + 8);   // 10
  EXPECT_EQ(matrixA[0][0], 0.0);
  EXPECT_EQ(matrixA[1][0], 1.0);
  EXPECT_EQ(matrixA[2][0], 2.0);
  EXPECT_EQ(matrixB[0][0], 0.0);
  EXPECT_EQ(matrixB[1][0], 2.0);
  EXPECT_EQ(matrixB[2][0], 4.0);

  return { matrixA, matrixB };
}

template<template<class> class MatrixPolicy>
void testMismatchedRowDimensions()
{
  MatrixPolicy<double> matrixA{ 3, 3, 1.0 };
  MatrixPolicy<double> matrixB{ 4, 3, 2.0 };  // Different number of rows!

  // Should throw when creating the function with mismatched row counts
  EXPECT_ANY_THROW(MatrixPolicy<double>::Function(
    [](auto&& mA, auto&& mB)
    {
      // This should throw when matrixA and matrixB have different row counts
      mA.ForEachRow([&](const double& a, const double& b, double& c)
        { c = a + b; },
        mA.GetConstColumnView(0),
        mB.GetConstColumnView(0),
        mA.GetColumnView(1));
    }, matrixA, matrixB));
}

template<template<class> class MatrixPolicy>
void testMismatchedColumnDimensions()
{
  MatrixPolicy<double> matrix{ 3, 4, 1.0 };

  // Create the function - this should succeed
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m)
    {
      // Try to access a column that doesn't exist
      m.ForEachRow([&](const double& a, double& b)
        { b = a * 2.0; },
        m.GetConstColumnView(0),
        m.GetColumnView(5));  // Column 5 doesn't exist in a 4-column matrix
    }, matrix);

  // Should throw when invoking the function because column 5 doesn't exist
  EXPECT_ANY_THROW(func(matrix));
}

template<template<class> class MatrixPolicy>
void testWrongMatrixDimensions()
{
  MatrixPolicy<double> matrix1{ 3, 4, 1.0 };
  MatrixPolicy<double> matrix2{ 3, 5, 2.0 };  // Different column count

  // Create a function that expects 4 columns
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m)
    {
      m.ForEachRow([&](const double& a, double& b)
        { b = a * 2.0; },
        m.GetConstColumnView(0),
        m.GetColumnView(3));  // Column 3 exists in 4-column matrix
    }, matrix1);

  func(matrix1);  // Should work fine
  EXPECT_NO_THROW(func(matrix1));

  // Should throw when applied to matrix with wrong dimensions
  EXPECT_ANY_THROW(func(matrix2));
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testMultipleTemporaries()
{
  MatrixPolicy<double> matrix{ 4, 5, 0.0 };

  // Initialize first two columns
  for (std::size_t i = 0; i < matrix.NumRows(); ++i)
  {
    matrix[i][0] = static_cast<double>(i + 1);
    matrix[i][1] = static_cast<double>((i + 1) * 10);
  }

  auto func = MatrixPolicy<double>::Function(
    [](auto&& m)
    {
      // Use TWO temporaries for intermediate calculations
      auto tmp1 = m.GetRowVariable();
      auto tmp2 = m.GetRowVariable();

      // tmp1 = col0 * col1
      m.ForEachRow([&](const double& a, const double& b, double& t)
        { t = a * b; },
        m.GetConstColumnView(0),
        m.GetConstColumnView(1),
        tmp1);

      // tmp2 = col0 + col1
      m.ForEachRow([&](const double& a, const double& b, double& t)
        { t = a + b; },
        m.GetConstColumnView(0),
        m.GetConstColumnView(1),
        tmp2);

      // col2 = tmp1 + tmp2 (product + sum)
      m.ForEachRow([&](const double& t1, const double& t2, double& c)
        { c = t1 + t2; },
        tmp1,
        tmp2,
        m.GetColumnView(2));

      // col3 = tmp1 - tmp2 (product - sum)
      m.ForEachRow([&](const double& t1, const double& t2, double& c)
        { c = t1 - t2; },
        tmp1,
        tmp2,
        m.GetColumnView(3));

      // col4 = tmp1 * tmp2
      m.ForEachRow([&](const double& t1, const double& t2, double& c)
        { c = t1 * t2; },
        tmp1,
        tmp2,
        m.GetColumnView(4));
    }, matrix);

  func(matrix);

  // Verify results
  // Row 0: col0=1, col1=10, product=10, sum=11
  EXPECT_EQ(matrix[0][2], 10.0 + 11.0);   // 21
  EXPECT_EQ(matrix[0][3], 10.0 - 11.0);   // -1
  EXPECT_EQ(matrix[0][4], 10.0 * 11.0);   // 110

  // Row 1: col0=2, col1=20, product=40, sum=22
  EXPECT_EQ(matrix[1][2], 40.0 + 22.0);   // 62
  EXPECT_EQ(matrix[1][3], 40.0 - 22.0);   // 18
  EXPECT_EQ(matrix[1][4], 40.0 * 22.0);   // 880

  // Row 3: col0=4, col1=40, product=160, sum=44
  EXPECT_EQ(matrix[3][2], 160.0 + 44.0);  // 204
  EXPECT_EQ(matrix[3][3], 160.0 - 44.0);  // 116
  EXPECT_EQ(matrix[3][4], 160.0 * 44.0);  // 7040

  return matrix;
}

template<template<class> class MatrixPolicy>
MatrixPolicy<double> testColumnViewReuse()
{
  MatrixPolicy<double> matrix{ 3, 4, 0.0 };

  for (std::size_t i = 0; i < matrix.NumRows(); ++i)
    matrix[i][0] = static_cast<double>(i + 1);

  auto func = MatrixPolicy<double>::Function(
    [](auto&& m)
    {
      // Create column views once
      auto col0 = m.GetConstColumnView(0);
      auto col1 = m.GetColumnView(1);
      auto col2 = m.GetColumnView(2);
      auto col3 = m.GetColumnView(3);

      // Reuse the same column views in multiple ForEachRow calls
      // col1 = col0 * 2
      m.ForEachRow([&](const double& a, double& b)
        { b = a * 2.0; },
        col0, col1);

      // col2 = col0 + col1 (reusing col0 and col1)
      m.ForEachRow([&](const double& a, const double& b, double& c)
        { c = a + b; },
        col0, col1, col2);

      // col3 = col2 * col1 (reusing col1 and col2)
      m.ForEachRow([&](const double& a, const double& b, double& c)
        { c = a * b; },
        col2, col1, col3);
    }, matrix);

  func(matrix);

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
MatrixPolicy<double> testFunctionReusability()
{
  // Create a function once
  MatrixPolicy<double> matrix1{ 2, 3, 1.0 };
  
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m)
    {
      auto tmp = m.GetRowVariable();
      m.ForEachRow([&](const double& a, const double& b, const double& c, double& t)
        { t = a + b + c; },
        m.GetConstColumnView(0),
        m.GetConstColumnView(1),
        m.GetConstColumnView(2),
        tmp);
      m.ForEachRow([&](double& c, const double& t)
        { c = 2.0 * t; },
        m.GetColumnView(2),
        tmp);
    }, matrix1);

  // Apply to first matrix
  for (std::size_t i = 0; i < matrix1.NumRows(); ++i)
    for (std::size_t j = 0; j < matrix1.NumColumns(); ++j)
      matrix1[i][j] = static_cast<double>(i + j);

  func(matrix1);
  EXPECT_EQ(matrix1[0][2], 2.0 * (0 + 1 + 2));  // 6
  EXPECT_EQ(matrix1[1][2], 2.0 * (1 + 2 + 3));  // 12

  // Apply to second matrix with same dimensions
  MatrixPolicy<double> matrix2{ 2, 3, 5.0 };
  func(matrix2);
  EXPECT_EQ(matrix2[0][2], 2.0 * (5 + 5 + 5));  // 30
  EXPECT_EQ(matrix2[1][2], 2.0 * (5 + 5 + 5));  // 30

  // Apply to third matrix with different values
  MatrixPolicy<double> matrix3{ 2, 3, 0.0 };
  for (std::size_t i = 0; i < matrix3.NumRows(); ++i)
    matrix3[i][0] = static_cast<double>(i * 10);

  func(matrix3);
  EXPECT_EQ(matrix3[0][2], 2.0 * (0 + 0 + 0));   // 0
  EXPECT_EQ(matrix3[1][2], 2.0 * (10 + 0 + 0));  // 20

  return matrix1;
}

template<template<class> class MatrixPolicy>
void testConstMatrixFunction()
{
  MatrixPolicy<double> matrix{ 3, 4, 0.0 };
  
  // Set initial values
  for (std::size_t i = 0; i < matrix.NumRows(); ++i)
    for (std::size_t j = 0; j < matrix.NumColumns(); ++j)
      matrix[i][j] = static_cast<double>(i * 10 + j);
  
  // Create a const reference
  const MatrixPolicy<double>& const_matrix = matrix;
  
  // Create a function that only reads from the matrix
  auto read_func = MatrixPolicy<double>::Function(
    [](auto&& m)
    {
      auto tmp = m.GetRowVariable();
      // Only use GetConstColumnView - should work with const matrices
      m.ForEachRow([&](const double& a, const double& b, double& t)
        { t = a + b; },
        m.GetConstColumnView(0),
        m.GetConstColumnView(1),
        tmp);
      
      // Verify we can read the values (no writes to m)
      double sum = 0.0;
      m.ForEachRow([&sum](const double& val)
        { sum += val; },
        m.GetConstColumnView(2));
    }, const_matrix);
  
  // Should work fine with const matrix
  EXPECT_NO_THROW(read_func(const_matrix));
  
  // Verify original matrix unchanged
  EXPECT_EQ(matrix[0][0], 0.0);
  EXPECT_EQ(matrix[1][2], 12.0);
}

template<template<class> class MatrixPolicy>
void testEmptyMatrixFunction()
{
  // Test with 0 rows
  MatrixPolicy<double> empty_rows{ 0, 3, 1.0 };
  
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m)
    {
      // This should never execute
      m.ForEachRow([&](double& val)
        { val = 99.0; },  // Would fail if executed
        m.GetColumnView(0));
    }, empty_rows);
  
  // Should not throw, just iterate 0 times
  EXPECT_NO_THROW(func(empty_rows));
  
  // Test with 0 columns (edge case)
  MatrixPolicy<double> empty_cols{ 3, 0, 1.0 };
  
  auto func2 = MatrixPolicy<double>::Function(
    [](auto&&)
    {
      // Cannot get any column views, so just return
    }, empty_cols);
  
  EXPECT_NO_THROW(func2(empty_cols));
}

