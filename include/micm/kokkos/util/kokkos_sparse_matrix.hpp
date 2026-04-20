// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/sparse_matrix.hpp>

#include <Kokkos_Core.hpp>
#include <vector>

namespace micm
{
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
    ViewType view_;
    HostViewType h_view_;

   public:
    KokkosSparseMatrix()
        : SparseMatrix<T, OrderingPolicy>()
    {
    }

    KokkosSparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder, bool indexing_only = false)
        : SparseMatrix<T, OrderingPolicy>(builder, indexing_only)
    {
      view_ = ViewType("sparse_matrix", this->data_.size());
      h_view_ = Kokkos::create_mirror_view(view_);
    }

    void CopyToDevice()
    {
      if (view_.extent(0) != this->data_.size())
      {
        view_ = ViewType("sparse_matrix", this->data_.size());
        h_view_ = Kokkos::create_mirror_view(view_);
      }
      for (std::size_t i = 0; i < this->data_.size(); ++i)
      {
        h_view_(i) = this->data_[i];
      }
      Kokkos::deep_copy(view_, h_view_);
    }

    void CopyToHost()
    {
      if (view_.extent(0) != 0)
      {
        Kokkos::deep_copy(h_view_, view_);
        for (std::size_t i = 0; i < this->data_.size(); ++i)
        {
          this->data_[i] = h_view_(i);
        }
      }
    }

    ViewType GetView() const
    {
      return view_;
    }

    void Fill(T val)
    {
      if (view_.extent(0) != this->data_.size())
      {
        view_ = ViewType("sparse_matrix", this->data_.size());
        h_view_ = Kokkos::create_mirror_view(view_);
      }
      Kokkos::deep_copy(view_, val);
    }
  };
}  // namespace micm
