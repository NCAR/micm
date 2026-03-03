// Tests of common matrix functions
#include <gtest/gtest.h>

#include <array>
#include <vector>

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

  // Use the function with a different matrix with the same number of columns, but different number of rows,
  // to test that it works with different sizes
  MatrixPolicy<double> matrix2{ 3, 3, -1.0 };
  func(matrix2);
  EXPECT_EQ(matrix2[0][2], 4.0 * (-1 + -1 + -1));   // -12
  EXPECT_EQ(matrix2[1][2], 4.0 * (-1 + -1 + -1));   // -12
  EXPECT_EQ(matrix2[2][2], 4.0 * (-1 + -1 + -1));   // -12
  EXPECT_EQ(matrix2[0][0], -1.0);
  EXPECT_EQ(matrix2[1][0], -1.0);
  EXPECT_EQ(matrix2[2][0], -1.0);
  EXPECT_EQ(matrix2[0][1], -1.0);
  EXPECT_EQ(matrix2[1][1], -1.0);
  EXPECT_EQ(matrix2[2][1], -1.0);

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
MatrixPolicy<double> testVectorInMatrixFunction()
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

  // Create a vector that we will use in the function
  std::vector<double> vec(matrix.NumRows());

  // Set some initial values in the vector
  vec[0] = 100.0;
  vec[1] = 200.0;
  vec[2] = 300.0;
  vec[3] = 400.0;
  vec[4] = 500.0;

  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      auto tmp = m.GetRowVariable();
      m.ForEachRow([&](const double& a, const double& b, const double& c, double& t)
        { t = a + b + c; },
        m.GetConstColumnView(0),
        m.GetConstColumnView(1),
        m.GetConstColumnView(2),
        tmp);
      m.ForEachRow([&](double& c, const double& d, const double& t)
        { c = d * t; },
        m.GetColumnView(2),
        v,
        tmp);
    }, matrix, vec); // pass matrix so the type and dimensions are known by the function

  func(matrix, vec); // apply the function to the matrix

  // Check results
  EXPECT_EQ(matrix[0][2], 100.0 * (-2 + 8 + 18));   // 2400
  EXPECT_EQ(matrix[1][2], 200.0 * (-1 + 9 + 19));   // 5400
  EXPECT_EQ(matrix[2][2], 300.0 * (0 + 10 + 20));   // 9000
  EXPECT_EQ(matrix[3][2], 400.0 * (1 + 11 + 21));   // 13200
  EXPECT_EQ(matrix[4][2], 500.0 * (2 + 12 + 22));   // 17000
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
std::tuple<MatrixPolicy<double>, MatrixPolicy<double>> testMultiMatrixDifferentRowsFromCreation()
{
  // Create function with 3-row matrices
  MatrixPolicy<double> matrixA_create{ 3, 2, 0.0 };
  MatrixPolicy<double> matrixB_create{ 3, 3, 0.0 };
  
  auto func = MatrixPolicy<double>::Function(
    [](auto&& mA, auto&& mB)
    {
      auto tmp = mA.GetRowVariable();
      mA.ForEachRow([&](const double& a, const double& b, double& t)
        { t = a + b; },
        mA.GetConstColumnView(0),
        mB.GetConstColumnView(2),
        tmp);
      mA.ForEachRow([&](const double& t, double& c)
        { c = t * 2.0; },
        tmp,
        mA.GetColumnView(1));
    }, matrixA_create, matrixB_create);
  
  // Now use with 5-row matrices (different from creation)
  MatrixPolicy<double> matrixA{ 5, 2, 0.0 };
  MatrixPolicy<double> matrixB{ 5, 3, 0.0 };
  
  for (std::size_t i = 0; i < 5; ++i)
  {
    matrixA[i][0] = static_cast<double>(i + 1);
    matrixB[i][2] = static_cast<double>(i * 10);
  }
  
  // Should work - column counts match, row counts match each other
  func(matrixA, matrixB);
  
  EXPECT_EQ(matrixA[0][1], (1.0 + 0.0) * 2.0);    // 2
  EXPECT_EQ(matrixA[1][1], (2.0 + 10.0) * 2.0);   // 24
  EXPECT_EQ(matrixA[2][1], (3.0 + 20.0) * 2.0);   // 46
  EXPECT_EQ(matrixA[3][1], (4.0 + 30.0) * 2.0);   // 68
  EXPECT_EQ(matrixA[4][1], (5.0 + 40.0) * 2.0);   // 90
  
  // Also test with 2-row matrices (fewer rows than creation)
  MatrixPolicy<double> matrixA2{ 2, 2, 0.0 };
  MatrixPolicy<double> matrixB2{ 2, 3, 0.0 };
  
  matrixA2[0][0] = 10.0;
  matrixA2[1][0] = 20.0;
  matrixB2[0][2] = 5.0;
  matrixB2[1][2] = 15.0;
  
  func(matrixA2, matrixB2);
  
  EXPECT_EQ(matrixA2[0][1], (10.0 + 5.0) * 2.0);   // 30
  EXPECT_EQ(matrixA2[1][1], (20.0 + 15.0) * 2.0);  // 70
  
  return { matrixA, matrixB };
}

/// @brief Test: Matrix + vector - function created with N rows, used with M rows
template<template<class> class MatrixPolicy>
std::tuple<MatrixPolicy<double>, std::vector<double>> testMatrixVectorDifferentRowsFromCreation()
{
  // Create function with 3-row matrix and vector
  MatrixPolicy<double> matrix_create{ 3, 3, 0.0 };
  std::vector<double> vec_create(3);
  
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      auto tmp = m.GetRowVariable();
      m.ForEachRow([&](const double& a, const double& b, double& t)
        { t = a + b; },
        m.GetConstColumnView(0),
        v,
        tmp);
      m.ForEachRow([&](const double& t, double& c)
        { c = t * 3.0; },
        tmp,
        m.GetColumnView(1));
    }, matrix_create, vec_create);
  
  // Now use with 5-row matrix and vector (different from creation)
  MatrixPolicy<double> matrix{ 5, 3, 0.0 };
  std::vector<double> vec(5);
  
  for (std::size_t i = 0; i < 5; ++i)
  {
    matrix[i][0] = static_cast<double>(i + 1);
    vec[i] = static_cast<double>(i * 10);
  }
  
  // Should work - columns match, row counts match each other
  func(matrix, vec);
  
  EXPECT_EQ(matrix[0][1], (1.0 + 0.0) * 3.0);    // 3
  EXPECT_EQ(matrix[1][1], (2.0 + 10.0) * 3.0);   // 36
  EXPECT_EQ(matrix[2][1], (3.0 + 20.0) * 3.0);   // 69
  EXPECT_EQ(matrix[3][1], (4.0 + 30.0) * 3.0);   // 102
  EXPECT_EQ(matrix[4][1], (5.0 + 40.0) * 3.0);   // 135
  
  return { matrix, vec };
}

/// @brief Test: Mismatched row counts at invocation time (should fail)
template<template<class> class MatrixPolicy>
void testMismatchedRowsAtInvocation()
{
  MatrixPolicy<double> matrix_create{ 3, 2, 0.0 };
  std::vector<double> vec_create(3);
  
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      m.ForEachRow([&](const double& a, double& b)
        { b = a * 2.0; },
        v,
        m.GetColumnView(0));
    }, matrix_create, vec_create);
  
  // Try to invoke with matrix (5 rows) and vector (3 rows) - should fail
  MatrixPolicy<double> matrix{ 5, 2, 0.0 };
  std::vector<double> vec(3);
  
  EXPECT_ANY_THROW(func(matrix, vec));
  
  // Try the other way - matrix (3 rows) and vector (5 rows) - should also fail
  MatrixPolicy<double> matrix2{ 3, 2, 0.0 };
  std::vector<double> vec2(5);
  
  EXPECT_ANY_THROW(func(matrix2, vec2));
}

/// @brief Test: Mismatched row counts between multiple matrices at invocation (should fail)
template<template<class> class MatrixPolicy>
void testMultipleMatricesMismatchedRowsAtInvocation()
{
  MatrixPolicy<double> matrixA_create{ 3, 2, 0.0 };
  MatrixPolicy<double> matrixB_create{ 3, 3, 0.0 };
  
  auto func = MatrixPolicy<double>::Function(
    [](auto&& mA, auto&& mB)
    {
      mA.ForEachRow([&](const double& a, double& b)
        { b = a * 2.0; },
        mB.GetConstColumnView(0),
        mA.GetColumnView(0));
    }, matrixA_create, matrixB_create);
  
  // Try to invoke with matrices having different row counts - should fail
  MatrixPolicy<double> matrixA{ 5, 2, 0.0 };
  MatrixPolicy<double> matrixB{ 3, 3, 0.0 };  // Different row count!
  
  EXPECT_ANY_THROW(func(matrixA, matrixB));
}

/// @brief Test: Wrong column count at invocation time (should fail)
template<template<class> class MatrixPolicy>
void testWrongColumnCountAtInvocation()
{
  // Create function with 3-column matrix
  MatrixPolicy<double> matrix_create{ 4, 3, 0.0 };
  
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
    }, matrix_create);
  
  // Try to invoke with wrong column count - should fail
  MatrixPolicy<double> matrix_wrong_cols{ 4, 4, 0.0 };  // 4 columns instead of 3
  EXPECT_ANY_THROW(func(matrix_wrong_cols));
  
  MatrixPolicy<double> matrix_wrong_cols2{ 4, 2, 0.0 };  // 2 columns instead of 3
  EXPECT_ANY_THROW(func(matrix_wrong_cols2));
  
  // Should work with different row count but same column count
  MatrixPolicy<double> matrix_ok{ 7, 3, 0.0 };  // 7 rows, 3 columns
  EXPECT_NO_THROW(func(matrix_ok));
}

template<template<class> class MatrixPolicy>
void testMismatchedRowDimensions()
{
  MatrixPolicy<double> matrixA{ 3, 3, 1.0 };
  MatrixPolicy<double> matrixB{ 4, 3, 2.0 };  // Different number of rows during creation

  // Should now SUCCEED when creating with different row counts
  // (as long as column counts match, which they do here)
  auto func = MatrixPolicy<double>::Function(
    [](auto&& mA, auto&& mB)
    {
      mA.ForEachRow([&](const double& a, const double& b, double& c)
        { c = a + b; },
        mA.GetConstColumnView(0),
        mB.GetConstColumnView(0),
        mA.GetColumnView(1));
    }, matrixA, matrixB);
  
  // Can use the function with matrices of the same row count
  MatrixPolicy<double> matrixC{ 5, 3, 0.0 };
  MatrixPolicy<double> matrixD{ 5, 3, 1.0 };
  
  EXPECT_NO_THROW(func(matrixC, matrixD));
  
  // But should throw if invoked with mismatched row counts
  MatrixPolicy<double> matrixE{ 3, 3, 0.0 };
  MatrixPolicy<double> matrixF{ 4, 3, 1.0 };  // Different row count!
  
  EXPECT_ANY_THROW(func(matrixE, matrixF));
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

  // Should throw when applied to matrix with wrong column count
  EXPECT_ANY_THROW(func(matrix2));
  
  // But should work with different row count as long as column count matches
  MatrixPolicy<double> matrix3{ 7, 4, 1.0 };  // 7 rows, 4 columns
  EXPECT_NO_THROW(func(matrix3));
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

// ============================================================================
// Vector Support Validation Tests
// ============================================================================

/// @brief Test: Vector with TOO FEW elements (creation succeeds, invocation fails)
template<template<class> class MatrixPolicy>
void testVectorTooSmall()
{
  MatrixPolicy<double> matrix{ 5, 3, 1.0 };
  std::vector<double> vec_too_small(3);  // Only 3 elements, but matrix has 5 rows
  
  // Should succeed at creation (row counts can differ at creation)
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      auto tmp = m.GetRowVariable();
      m.ForEachRow([&](const double& a, const double& b, double& t)
        { t = a + b; },
        m.GetConstColumnView(0),
        v,
        tmp);
    }, matrix, vec_too_small);
  
  // Should throw at invocation when vector size doesn't match matrix row count
  EXPECT_ANY_THROW(func(matrix, vec_too_small));
}

/// @brief Test: Vector with TOO MANY elements (creation succeeds, invocation fails)
template<template<class> class MatrixPolicy>
void testVectorTooLarge()
{
  MatrixPolicy<double> matrix{ 5, 3, 1.0 };
  std::vector<double> vec_too_large(10);  // 10 elements, but matrix has 5 rows
  
  // Should succeed at creation (row counts can differ at creation)
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      auto tmp = m.GetRowVariable();
      m.ForEachRow([&](const double& a, const double& b, double& t)
        { t = a + b; },
        m.GetConstColumnView(0),
        v,
        tmp);
    }, matrix, vec_too_large);
  
  // Should throw at invocation when vector size doesn't match matrix row count
  EXPECT_ANY_THROW(func(matrix, vec_too_large));
}

/// @brief Test: Empty vector with non-empty matrix (creation succeeds, invocation fails)
template<template<class> class MatrixPolicy>
void testEmptyVectorNonEmptyMatrix()
{
  MatrixPolicy<double> matrix{ 5, 3, 1.0 };
  std::vector<double> empty_vec;  // Empty
  
  // Should succeed at creation
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      m.ForEachRow([&](const double& a, double& b)
        { b = a; },
        v,
        m.GetColumnView(0));
    }, matrix, empty_vec);
  
  // Should throw at invocation when vector size doesn't match
  EXPECT_ANY_THROW(func(matrix, empty_vec));
}

/// @brief Test: Non-empty vector with empty matrix (creation succeeds, invocation fails)
template<template<class> class MatrixPolicy>
void testNonEmptyVectorEmptyMatrix()
{
  MatrixPolicy<double> matrix{ 0, 3, 1.0 };  // 0 rows
  std::vector<double> vec(5);
  
  // Should succeed at creation
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      m.ForEachRow([&](const double& a, double& b)
        { b = a; },
        v,
        m.GetColumnView(0));
    }, matrix, vec);
  
  // Should throw at invocation when vector size doesn't match matrix row count
  EXPECT_ANY_THROW(func(matrix, vec));
}

/// @brief Test: Empty vector with empty matrix (should work - no iterations)
template<template<class> class MatrixPolicy>
void testEmptyVectorEmptyMatrix()
{
  MatrixPolicy<double> matrix{ 0, 3, 1.0 };  // 0 rows
  std::vector<double> empty_vec;  // Empty
  
  // Should succeed - both are empty, ForEachRow won't iterate
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      m.ForEachRow([&](const double& a, double& b)
        { b = a; },
        v,
        m.GetColumnView(0));
    }, matrix, empty_vec);
  
  EXPECT_NO_THROW(func(matrix, empty_vec));
}

/// @brief Test: Multiple vectors with DIFFERENT sizes (creation succeeds, invocation fails)
template<template<class> class MatrixPolicy>
void testMultipleVectorsDifferentSizes()
{
  MatrixPolicy<double> matrix{ 5, 3, 1.0 };
  std::vector<double> vec1(5);  // Size 5
  std::vector<double> vec2(3);  // Size 3 - different!
  
  // Should succeed at creation (different row counts allowed at creation)
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v1, auto&& v2)
    {
      m.ForEachRow([&](const double& a, const double& b, double& c)
        { c = a + b; },
        v1,
        v2,
        m.GetColumnView(0));
    }, matrix, vec1, vec2);
  
  // Should throw at invocation because vectors have different sizes
  EXPECT_ANY_THROW(func(matrix, vec1, vec2));
}

/// @brief Test: Multiple vectors with SAME correct size (should work)
template<template<class> class MatrixPolicy>
MatrixPolicy<double> testMultipleVectorsSameSize()
{
  MatrixPolicy<double> matrix{ 5, 3, 0.0 };
  std::vector<double> vec1(5);
  std::vector<double> vec2(5);
  
  // Initialize vectors
  for (std::size_t i = 0; i < 5; ++i)
  {
    vec1[i] = static_cast<double>(i + 1);
    vec2[i] = static_cast<double>((i + 1) * 10);
  }
  
  // Should succeed
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v1, auto&& v2)
    {
      // col0 = v1 + v2
      m.ForEachRow([&](const double& a, const double& b, double& c)
        { c = a + b; },
        v1,
        v2,
        m.GetColumnView(0));
    }, matrix, vec1, vec2);
  
  func(matrix, vec1, vec2);
  
  // Verify results
  EXPECT_EQ(matrix[0][0], 1.0 + 10.0);    // 11
  EXPECT_EQ(matrix[1][0], 2.0 + 20.0);    // 22
  EXPECT_EQ(matrix[2][0], 3.0 + 30.0);    // 33
  EXPECT_EQ(matrix[3][0], 4.0 + 40.0);    // 44
  EXPECT_EQ(matrix[4][0], 5.0 + 50.0);    // 55
  
  return matrix;
}

/// @brief Test: Multiple matrices + vector - vector size must match all matrices
template<template<class> class MatrixPolicy>
std::tuple<MatrixPolicy<double>, MatrixPolicy<double>> testMultipleMatricesOneVector()
{
  MatrixPolicy<double> matrixA{ 4, 2, 0.0 };
  MatrixPolicy<double> matrixB{ 4, 3, 0.0 };
  std::vector<double> vec(4);
  
  for (std::size_t i = 0; i < 4; ++i)
  {
    vec[i] = static_cast<double>(i * 2);
    matrixA[i][0] = static_cast<double>(i + 1);
    matrixB[i][0] = static_cast<double>(i * 10);
  }
  
  // Should succeed - both matrices have 4 rows, vector has 4 elements
  auto func = MatrixPolicy<double>::Function(
    [](auto&& mA, auto&& mB, auto&& v)
    {
      // matrixA col1 = matrixA col0 + vector
      mA.ForEachRow([&](const double& a, const double& b, double& c)
        { c = a + b; },
        mA.GetConstColumnView(0),
        v,
        mA.GetColumnView(1));
      
      // matrixB col1 = matrixB col0 + vector
      mB.ForEachRow([&](const double& a, const double& b, double& c)
        { c = a + b; },
        mB.GetConstColumnView(0),
        v,
        mB.GetColumnView(1));
    }, matrixA, matrixB, vec);
  
  func(matrixA, matrixB, vec);
  
  // Verify results
  EXPECT_EQ(matrixA[0][1], 1.0 + 0.0);   // 1
  EXPECT_EQ(matrixA[1][1], 2.0 + 2.0);   // 4
  EXPECT_EQ(matrixA[2][1], 3.0 + 4.0);   // 7
  EXPECT_EQ(matrixA[3][1], 4.0 + 6.0);   // 10
  
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
void testMultipleMatricesDifferentRowsVector()
{
  MatrixPolicy<double> matrixA{ 4, 2, 0.0 };
  MatrixPolicy<double> matrixB{ 5, 3, 0.0 };  // Different row count!
  std::vector<double> vec(4);
  
  // Should succeed at creation (different row counts allowed at creation)
  auto func = MatrixPolicy<double>::Function(
    [](auto&& mA, auto&& mB, auto&& v)
    {
      mA.ForEachRow([&](const double& a, double& b)
        { b = a; },
        v,
        mA.GetColumnView(0));
    }, matrixA, matrixB, vec);
  
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
void testVectorSizeMatchesOneMatrixOnly()
{
  MatrixPolicy<double> matrixA{ 5, 2, 0.0 };
  MatrixPolicy<double> matrixB{ 5, 3, 0.0 };
  std::vector<double> vec(4);  // Wrong size for both matrices (they have 5 rows)
  
  // Should succeed at creation (different row counts allowed at creation)
  auto func = MatrixPolicy<double>::Function(
    [](auto&& mA, auto&& mB, auto&& v)
    {
      mA.ForEachRow([&](const double& a, double& b)
        { b = a; },
        v,
        mA.GetColumnView(0));
    }, matrixA, matrixB, vec);
  
  // Should throw at invocation because vector size doesn't match matrix row counts
  EXPECT_ANY_THROW(func(matrixA, matrixB, vec));
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

/// @brief Test: Const vector (read-only access)
template<template<class> class MatrixPolicy>
MatrixPolicy<double> testConstVector()
{
  MatrixPolicy<double> matrix{ 3, 2, 0.0 };
  std::vector<double> vec_data = { 10.0, 20.0, 30.0 };
  const std::vector<double>& const_vec = vec_data;
  
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      // Read from const vector, write to matrix
      m.ForEachRow([&](const double& a, double& b)
        { b = a * 2.0; },
        v,  // const vector
        m.GetColumnView(0));
    }, matrix, const_vec);
  
  func(matrix, const_vec);
  
  EXPECT_EQ(matrix[0][0], 20.0);
  EXPECT_EQ(matrix[1][0], 40.0);
  EXPECT_EQ(matrix[2][0], 60.0);
  
  return matrix;
}

/// @brief Test: Non-const vector that gets modified
template<template<class> class MatrixPolicy>
std::tuple<MatrixPolicy<double>, std::vector<double>> testMutableVector()
{
  MatrixPolicy<double> matrix{ 3, 2, 0.0 };
  std::vector<double> vec = { 5.0, 10.0, 15.0 };
  
  for (std::size_t i = 0; i < 3; ++i)
    matrix[i][0] = static_cast<double>(i + 1);
  
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      // Write to vector from matrix
      m.ForEachRow([&](const double& a, double& b)
        { b = a * 3.0; },
        m.GetConstColumnView(0),
        v);  // non-const vector
    }, matrix, vec);
  
  func(matrix, vec);
  
  // Vector should be modified
  EXPECT_EQ(vec[0], 3.0);
  EXPECT_EQ(vec[1], 6.0);
  EXPECT_EQ(vec[2], 9.0);
  
  return { matrix, vec };
}

/// @brief Test: Function reusability with different vectors
template<template<class> class MatrixPolicy>
MatrixPolicy<double> testFunctionReusabilityWithVectors()
{
  MatrixPolicy<double> matrix{ 3, 2, 0.0 };
  std::vector<double> vec1 = { 1.0, 2.0, 3.0 };
  
  // Create function once
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      m.ForEachRow([&](const double& a, double& b)
        { b = a * 10.0; },
        v,
        m.GetColumnView(0));
    }, matrix, vec1);
  
  // Apply with first vector
  func(matrix, vec1);
  EXPECT_EQ(matrix[0][0], 10.0);
  EXPECT_EQ(matrix[1][0], 20.0);
  EXPECT_EQ(matrix[2][0], 30.0);
  
  // Apply with second vector (same size)
  std::vector<double> vec2 = { 5.0, 6.0, 7.0 };
  func(matrix, vec2);
  EXPECT_EQ(matrix[0][0], 50.0);
  EXPECT_EQ(matrix[1][0], 60.0);
  EXPECT_EQ(matrix[2][0], 70.0);
  
  return matrix;
}

/// @brief Test: Applying function with vector of wrong size at invocation time
template<template<class> class MatrixPolicy>
void testFunctionInvocationWithWrongSizedVector()
{
  MatrixPolicy<double> matrix{ 3, 2, 0.0 };
  std::vector<double> vec_correct(3);
  std::vector<double> vec_wrong(5);
  
  // Create function with correct-sized vector
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      m.ForEachRow([&](const double& a, double& b)
        { b = a; },
        v,
        m.GetColumnView(0));
    }, matrix, vec_correct);
  
  // Should work with correct size
  EXPECT_NO_THROW(func(matrix, vec_correct));
  
  // Should throw when invoked with wrong-sized vector
  EXPECT_ANY_THROW(func(matrix, vec_wrong));
}

/// @brief Test: Array instead of vector (should work with any operator[] type)
template<template<class> class MatrixPolicy>
MatrixPolicy<double> testArraySupport()
{
  MatrixPolicy<double> matrix{ 4, 2, 0.0 };
  std::array<double, 4> arr = { 100.0, 200.0, 300.0, 400.0 };
  
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& a)
    {
      m.ForEachRow([&](const double& val, double& result)
        { result = val / 2.0; },
        a,
        m.GetColumnView(0));
    }, matrix, arr);
  
  func(matrix, arr);
  
  EXPECT_EQ(matrix[0][0], 50.0);
  EXPECT_EQ(matrix[1][0], 100.0);
  EXPECT_EQ(matrix[2][0], 150.0);
  EXPECT_EQ(matrix[3][0], 200.0);
  
  return matrix;
}

/// @brief Test: Mixed - vector, column view, and row variable together
template<template<class> class MatrixPolicy>
MatrixPolicy<double> testMixedVectorColumnViewRowVariable()
{
  MatrixPolicy<double> matrix{ 4, 3, 0.0 };
  std::vector<double> vec(4);
  
  for (std::size_t i = 0; i < 4; ++i)
  {
    matrix[i][0] = static_cast<double>(i + 1);
    matrix[i][1] = static_cast<double>((i + 1) * 10);
    vec[i] = static_cast<double>((i + 1) * 100);
  }
  
  auto func = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v)
    {
      auto tmp = m.GetRowVariable();
      
      // tmp = col0 + col1
      m.ForEachRow([&](const double& a, const double& b, double& t)
        { t = a + b; },
        m.GetConstColumnView(0),
        m.GetConstColumnView(1),
        tmp);
      
      // col2 = tmp + vector
      m.ForEachRow([&](const double& t, const double& v_elem, double& result)
        { result = t + v_elem; },
        tmp,
        v,
        m.GetColumnView(2));
    }, matrix, vec);
  
  func(matrix, vec);
  
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

/// @brief Test: Different integer types (std::vector<int>)
template<template<class> class MatrixPolicy>
MatrixPolicy<int> testIntegerVector()
{
  MatrixPolicy<int> matrix{ 3, 2, 0 };
  std::vector<int> vec = { 10, 20, 30 };
  
  auto func = MatrixPolicy<int>::Function(
    [](auto&& m, auto&& v)
    {
      m.ForEachRow([&](const int& a, int& b)
        { b = a * 2; },
        v,
        m.GetColumnView(0));
    }, matrix, vec);
  
  func(matrix, vec);
  
  EXPECT_EQ(matrix[0][0], 20);
  EXPECT_EQ(matrix[1][0], 40);
  EXPECT_EQ(matrix[2][0], 60);
  
  return matrix;
}

template<template<class> class MatrixPolicy>
void testFunctionWithConstSignature()
{
  MatrixPolicy<double> matrix{ 3, 2, 0.0 };
  std::vector<double> vec(3);
  
  // Create function
  auto func_auto = MatrixPolicy<double>::Function(
    [](auto&& m, auto&& v) {
      m.ForEachRow([&](const double& a, double& b) { b = a * 2.0; },
                   v, m.GetColumnView(0));
    }, matrix, vec);
  
  // Try to wrap in std::function with const signature
  std::function<void(MatrixPolicy<double>&, const std::vector<double>&)> func_std = func_auto;
  
  func_std(matrix, vec);
}

