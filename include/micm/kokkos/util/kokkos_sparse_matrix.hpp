// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/sparse_matrix.hpp>

#include <Kokkos_Core.hpp>
#include <vector>

namespace micm
{
  /// @brief Provides a Kokkos implementation to the SparseMatrix functionality.
  ///
  /// Inherits from SparseMatrix (the MICM host-side data layout) and maintains
  /// a Kokkos::View as a device-side mirror. The caller must explicitly call
  /// CopyToDevice() / CopyToHost() to synchronize, matching the CUDA matrix pattern.
  template<class T = double, class OrderingPolicy = SparseMatrixStandardOrdering>
  class KokkosSparseMatrix : public SparseMatrix<T, OrderingPolicy>
  {
   public:
    static constexpr std::size_t GroupVectorSize()
    {
      return OrderingPolicy::GroupVectorSize();
    }
    using value_type = T;
    using ViewType = Kokkos::View<T*>;
    using HostViewType = typename ViewType::host_mirror_type;

   private:
    /// Device-side (or unified) view — the Kokkos mirror of MICM's data_
    ViewType view_;

   public:
    KokkosSparseMatrix()
        : SparseMatrix<T, OrderingPolicy>()
    {
    }

    KokkosSparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder, bool indexing_only = false)
        : SparseMatrix<T, OrderingPolicy>(builder, indexing_only),
          view_("sparse_matrix", SparseMatrix<T, OrderingPolicy>(builder, indexing_only).AsVector().size())
    {
    }

    /// @brief Copy host data (MICM's data_) to the device view
    void CopyToDevice()
    {
      if (view_.extent(0) != this->data_.size())
      {
        view_ = ViewType("sparse_matrix", this->data_.size());
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
        view_ = ViewType("sparse_matrix", this->data_.size());
      }
      Kokkos::deep_copy(view_, val);
    }
  };
}  // namespace micm
