// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Linear combinations of the Rosenbrock stage vectors in one pass over the data.
//
// The Rosenbrock step builds each stage input, each stage right-hand side, the new state, and the error estimate as
// sums of earlier stage vectors. A Copy() and one Axpy() for each term read and write the output once for each term.
// These functions read each input once and write each output once. The terms are added in the same order as the
// Axpy() calls, so the results are the same.
#pragma once

#include <micm/util/matrix.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <array>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace micm
{
  /// @brief True when the dense matrix keeps its data in host memory, so that the solver can use the one-pass
  ///        functions. The Kokkos and CUDA matrices derive from VectorMatrix, but they keep their data on the device.
  template<class DenseMatrixPolicy>
  struct HasHostLinearCombination : std::false_type
  {
  };
  template<class T, Index L>
  struct HasHostLinearCombination<VectorMatrix<T, L>> : std::true_type
  {
  };
  template<class T>
  struct HasHostLinearCombination<Matrix<T>> : std::true_type
  {
  };
  template<class DenseMatrixPolicy>
  constexpr bool HAS_HOST_LINEAR_COMBINATION_V = HasHostLinearCombination<DenseMatrixPolicy>::value;

  namespace detail
  {
    /// @brief The maximum number of terms, equal to the maximum number of Rosenbrock stages
    constexpr Index MAX_LINEAR_COMBINATION_TERMS = 6;

    template<Index N, class T>
    void LinearCombination(
        T* out,
        const T* base,
        const std::array<Real, N>& coefficients,
        const std::array<const T*, N>& x,
        const Index size)
    {
      for (Index i = 0; i < size; ++i)
      {
        T sum = base[i];
        for (Index j = 0; j < N; ++j)
        {
          sum += coefficients[j] * x[j][i];
        }
        out[i] = sum;
      }
    }

    template<Index N, class T>
    void LinearCombination2(
        T* out_a,
        const T* base_a,
        const Real* coefficients_a,
        T* out_b,
        const Real* coefficients_b,
        const std::array<const T*, N>& x,
        const Index size)
    {
      std::array<Real, N> c_a{};
      std::array<Real, N> c_b{};
      for (Index j = 0; j < N; ++j)
      {
        c_a[j] = coefficients_a[j];
        c_b[j] = coefficients_b[j];
      }
      for (Index i = 0; i < size; ++i)
      {
        T sum_a = base_a[i];
        T sum_b{};
        for (Index j = 0; j < N; ++j)
        {
          const T xj = x[j][i];
          sum_a += c_a[j] * xj;
          sum_b += c_b[j] * xj;
        }
        out_a[i] = sum_a;
        out_b[i] = sum_b;
      }
    }

    template<Index N, class Matrix>
    std::array<const typename Matrix::value_type*, N> Pointers(const std::vector<Matrix>& x)
    {
      std::array<const typename Matrix::value_type*, N> pointers{};
      for (Index j = 0; j < N; ++j)
      {
        pointers[j] = x[j].AsVector().data();
      }
      return pointers;
    }

    /// @brief Calls f with the number of terms as a compile-time constant, so that the compiler unrolls the term loop
    template<class F>
    void DispatchTermCount(const Index count, F&& f)
    {
      switch (count)
      {
        case 1: f(std::integral_constant<Index, 1>{}); break;
        case 2: f(std::integral_constant<Index, 2>{}); break;
        case 3: f(std::integral_constant<Index, 3>{}); break;
        case 4: f(std::integral_constant<Index, 4>{}); break;
        case 5: f(std::integral_constant<Index, 5>{}); break;
        case 6: f(std::integral_constant<Index, 6>{}); break;
        default: throw std::invalid_argument("A linear combination has from 1 to 6 terms");
      }
    }
  }  // namespace detail

  /// @brief Calculates out = base + sum of coefficients[j] * x[j] for j < count, in one pass
  /// @param out The output matrix. It can be the same matrix as base.
  /// @param base The matrix to add the terms to
  /// @param coefficients The coefficients of the terms
  /// @param x The matrices of the terms. Only the first count matrices are used.
  /// @param count The number of terms
  template<class Matrix>
  void LinearCombination(
      Matrix& out,
      const Matrix& base,
      const Real* coefficients,
      const std::vector<Matrix>& x,
      const Index count)
  {
    detail::DispatchTermCount(
        count,
        [&](auto n)
        {
          constexpr Index N = decltype(n)::value;
          std::array<Real, N> c{};
          for (Index j = 0; j < N; ++j)
          {
            c[j] = coefficients[j];
          }
          detail::LinearCombination<N>(
              out.AsVector().data(), base.AsVector().data(), c, detail::Pointers<N>(x), out.AsVector().size());
        });
  }

  /// @brief Calculates out_a = base_a + sum of coefficients_a[j] * x[j] and out_b = sum of coefficients_b[j] * x[j]
  ///        for j < count, in one pass
  template<class Matrix>
  void LinearCombination2(
      Matrix& out_a,
      const Matrix& base_a,
      const Real* coefficients_a,
      Matrix& out_b,
      const Real* coefficients_b,
      const std::vector<Matrix>& x,
      const Index count)
  {
    detail::DispatchTermCount(
        count,
        [&](auto n)
        {
          constexpr Index N = decltype(n)::value;
          detail::LinearCombination2<N>(
              out_a.AsVector().data(),
              base_a.AsVector().data(),
              coefficients_a,
              out_b.AsVector().data(),
              coefficients_b,
              detail::Pointers<N>(x),
              out_a.AsVector().size());
        });
  }
}  // namespace micm
