// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include <micm/util/view_category.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/vector_matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <vector>
#include <array>

using namespace micm;

// Test that ViewCategory SFINAE works correctly with the fallback specialization
TEST(ViewCategorySFINAE, VectorTypes)
{
  // std::vector should have ViewCategory_t = void (fallback)
  static_assert(std::same_as<ViewCategory_t<std::vector<double>>, void>);
  static_assert(std::same_as<ViewCategory_t<std::vector<int>>, void>);
  
  // std::array should also have ViewCategory_t = void (fallback)
  static_assert(std::same_as<ViewCategory_t<std::array<double, 10>>, void>);
  
  // VectorLike concept should accept vector types
  static_assert(VectorLike<std::vector<double>>);
  static_assert(VectorLike<std::vector<int>>);
  static_assert(VectorLike<std::array<double, 10>>);
}

TEST(ViewCategorySFINAE, MatrixViews)
{
  // Matrix column views should have DenseMatrixColumnViewTag
  using MatrixColumnView = Matrix<double>::ColumnView;
  static_assert(std::same_as<ViewCategory_t<MatrixColumnView>, DenseMatrixColumnViewTag>);
  static_assert(DenseMatrixColumnView<MatrixColumnView>);
  static_assert(!VectorLike<MatrixColumnView>);
  
  // VectorMatrix column views should have DenseMatrixColumnViewTag
  using VectorMatrixColumnView = VectorMatrix<double, 4>::ColumnView;
  static_assert(std::same_as<ViewCategory_t<VectorMatrixColumnView>, DenseMatrixColumnViewTag>);
  static_assert(DenseMatrixColumnView<VectorMatrixColumnView>);
  static_assert(!VectorLike<VectorMatrixColumnView>);
}

TEST(ViewCategorySFINAE, SparseMatrixViews)
{
  // Sparse matrix block views should have SparseMatrixBlockViewTag
  using SparseBlockView = SparseMatrix<double>::BlockView;
  static_assert(std::same_as<ViewCategory_t<SparseBlockView>, SparseMatrixBlockViewTag>);
  static_assert(SparseMatrixBlockView<SparseBlockView>);
  static_assert(!VectorLike<SparseBlockView>);
  
  using SparseConstBlockView = SparseMatrix<double>::ConstBlockView;
  static_assert(std::same_as<ViewCategory_t<SparseConstBlockView>, SparseMatrixBlockViewTag>);
  static_assert(SparseMatrixBlockView<SparseConstBlockView>);
  static_assert(!VectorLike<SparseConstBlockView>);
}

TEST(ViewCategorySFINAE, MatrixTypesExcluded)
{
  // Matrix types themselves should be excluded from VectorLike
  static_assert(!VectorLike<Matrix<double>>);
  static_assert(!VectorLike<VectorMatrix<double, 4>>);
  static_assert(!VectorLike<SparseMatrix<double>>);
}

TEST(ViewCategorySFINAE, RuntimeVectorUse)
{
  // Verify we can actually use std::vector in a VectorLike context
  std::vector<double> vec = {1.0, 2.0, 3.0, 4.0, 5.0};
  
  auto process_vector = []<VectorLike V>(V& v) -> double {
    double sum = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i) {
      sum += v[i];
    }
    return sum;
  };
  
  double sum = process_vector(vec);
  EXPECT_DOUBLE_EQ(sum, 15.0);
}
