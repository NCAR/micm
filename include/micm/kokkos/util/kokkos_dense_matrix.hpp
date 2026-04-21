// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/vector_matrix.hpp>

#include <Kokkos_Core.hpp>
#include <type_traits>
#include <vector>

namespace micm
{
  /// Concept for Kokkos Matrix
  template<typename MatrixType>
  concept KokkosMatrix = requires(MatrixType t) {
    { t.CopyToDevice() } -> std::same_as<void>;
    { t.CopyToHost() } -> std::same_as<void>;
  };

  /// @brief Provides a Kokkos implementation to the VectorMatrix functionality.
  ///
  /// Inherits from VectorMatrix (the MICM host-side data layout) and maintains
  /// a Kokkos::View as a device-side mirror. The caller must explicitly call
  /// CopyToDevice() / CopyToHost() to synchronize, matching the CUDA matrix pattern.
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
    using HostViewType = typename ViewType::host_mirror_type;

   private:
    /// Device-side (or unified) view — the Kokkos mirror of MICM's data_
    ViewType view_;

   public:
    KokkosDenseMatrix()
        : VectorMatrix<T, L>()
    {
    }

    KokkosDenseMatrix(std::size_t x_dim, std::size_t y_dim)
        : VectorMatrix<T, L>(x_dim, y_dim),
          view_("dense_matrix", VectorMatrix<T, L>(x_dim, y_dim).AsVector().size())
    {
    }

    KokkosDenseMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value)
        : VectorMatrix<T, L>(x_dim, y_dim, initial_value),
          view_("dense_matrix", VectorMatrix<T, L>(x_dim, y_dim).AsVector().size())
    {
      Kokkos::deep_copy(view_, initial_value);
    }

    /// @brief Copy host data (MICM's data_) to the device view
    void CopyToDevice()
    {
      if (view_.extent(0) != this->data_.size())
      {
        view_ = ViewType("dense_matrix", this->data_.size());
      }
      auto h_view = Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
          this->data_.data(), this->data_.size());
      Kokkos::deep_copy(view_, h_view);
    }

    /// @brief Copy device view data back to host (MICM's data_)
    void CopyToHost()
    {
      if (view_.extent(0) != 0)
      {
        auto h_view = Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            this->data_.data(), this->data_.size());
        Kokkos::deep_copy(h_view, view_);
      }
    }

    ViewType GetView() const
    {
      return view_;
    }

    /// @brief Set every element on the device to a given value
    void Fill(T val)
    {
      if (view_.extent(0) != this->data_.size())
      {
        view_ = ViewType("dense_matrix", this->data_.size());
      }
      Kokkos::deep_copy(view_, val);
    }
  };
}  // namespace micm
