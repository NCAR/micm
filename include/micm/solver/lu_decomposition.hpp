// Copyright (C) 2023-2025 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "lu_decomposition_doolittle.hpp"
#include "lu_decomposition_doolittle_in_place.hpp"
#include "lu_decomposition_mozart.hpp"
#include "lu_decomposition_mozart_in_place.hpp"

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

namespace micm
{

  /// @brief Concept for in-place LU decomposition algorithms
  template<class T, class SparseMatrixPolicy>
  concept LuDecompositionInPlaceConcept = requires(T t) {
    { t.GetLUMatrix(SparseMatrixPolicy{}, 0.0) };
    { t.Decompose(std::declval<SparseMatrixPolicy&>()) };
  };
  static_assert(
      LuDecompositionInPlaceConcept<LuDecompositionMozartInPlace, StandardSparseMatrix>,
      "LuDecompositionMozartInPlace does not meet the LuDecompositionInPlaceConcept requirements");
  static_assert(
      LuDecompositionInPlaceConcept<
          LuDecompositionMozartInPlace,
          SparseMatrix<double, SparseMatrixVectorOrderingCompressedSparseRow<1>>>,
      "LuDecompositionMozartInPlace for vector matrices does not meet the LuDecompositionInPlaceConcept requirements");
  static_assert(
      LuDecompositionInPlaceConcept<LuDecompositionDoolittleInPlace, StandardSparseMatrix>,
      "LuDecompositionDoolittleInPlace does not meet the LuDecompositionInPlaceConcept requirements");
  static_assert(
      LuDecompositionInPlaceConcept<
          LuDecompositionDoolittleInPlace,
          SparseMatrix<double, SparseMatrixVectorOrderingCompressedSparseRow<1>>>,
      "LuDecompositionDoolittleInPlace for vector matrices does not meet the LuDecompositionInPlaceConcept requirements");

  /// @brief Alias for the default LU decomposition algorithm
  using LuDecomposition = LuDecompositionDoolittle;

  /// @brief Alias for the default in-place LU decomposition algorithm
  using LuDecompositionInPlace = LuDecompositionMozartInPlace;
}  // namespace micm