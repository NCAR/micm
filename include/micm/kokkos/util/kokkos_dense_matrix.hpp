// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/kokkos/util/kokkos_util.hpp>
#include <micm/util/vector_matrix.hpp>

#include <Kokkos_Core.hpp>
#include <type_traits>
#include <vector>

namespace micm
{
  /**
   * @brief Provides a Kokkos implementation to the VectorMatrix functionality.
   */

  /// Concept for Kokkos Matrix
  template<typename MatrixType>
  concept KokkosMatrix = requires(MatrixType t) {
    { t.CopyToDevice() } -> std::same_as<void>;
    { t.CopyToHost() } -> std::same_as<void>;
  };

  template<class T, std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  class KokkosDenseMatrix : public VectorMatrix<T, L>
  {
   public:
    static constexpr std::size_t GroupVectorSize()
    {
      return L;
    }
    using value_type = T;
    using ViewType = Kokkos::View<T*>;
    using HostViewType = Kokkos::View<T*, Kokkos::HostSpace>;

   private:
    ViewType d_view_;
    HostViewType h_view_;

   public:
    KokkosDenseMatrix()
        : VectorMatrix<T, L>()
    {
      micm::kokkos::Initialize();
    }

    KokkosDenseMatrix(std::size_t x_dim, std::size_t y_dim)
        : VectorMatrix<T, L>(x_dim, y_dim)
    {
      micm::kokkos::Initialize();
    }

    KokkosDenseMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value)
        : VectorMatrix<T, L>(x_dim, y_dim, initial_value)
    {
      micm::kokkos::Initialize();
      Fill(initial_value);
    }

    void CopyToDevice()
    {
      if (d_view_.extent(0) != this->data_.size())
      {
        d_view_ = ViewType("dense_matrix", this->data_.size());
        h_view_ = Kokkos::create_mirror_view(d_view_);
      }
      for (std::size_t i = 0; i < this->data_.size(); ++i)
      {
        h_view_(i) = this->data_[i];
      }
      Kokkos::deep_copy(d_view_, h_view_);
    }

    void CopyToHost()
    {
      if (d_view_.extent(0) != 0)
      {
        Kokkos::deep_copy(h_view_, d_view_);
        for (std::size_t i = 0; i < this->data_.size(); ++i)
        {
          this->data_[i] = h_view_(i);
        }
      }
    }

    ViewType GetView() const
    {
      return d_view_;
    }

    void Fill(T val)
    {
      if (d_view_.extent(0) != this->data_.size())
      {
        d_view_ = ViewType("dense_matrix", this->data_.size());
        h_view_ = Kokkos::create_mirror_view(d_view_);
      }
      Kokkos::deep_copy(d_view_, val);
    }
  };
}  // namespace micm
