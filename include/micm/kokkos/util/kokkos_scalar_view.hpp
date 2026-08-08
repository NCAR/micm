// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

#include <Kokkos_Core.hpp>

namespace micm
{
  /// @brief A scalar view class for use in Matrix::Function lambdas
  template<class T>
  class KokkosScalarView
  {
    mutable T data_; // host data
    mutable Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> host_view_;
    mutable Kokkos::View<T*> device_view_{ "scalar", 1 }; // device data

   public:
    using value_type = T;
    struct View;
    struct ConstView;
    using ViewType = View;
    using ConstViewType = ConstView;

    struct ConstDeviceView
    {
        Kokkos::View<const T*> view_;

        KOKKOS_INLINE_FUNCTION operator T() const { return view_(0); }

        KOKKOS_INLINE_FUNCTION const T* data() const { return view_.data(); }
        
        KOKKOS_INLINE_FUNCTION ConstView GetView() const { return *this; }
    };

    struct DeviceView
    {
        Kokkos::View<T*> view_;

        KOKKOS_INLINE_FUNCTION operator T() const { return view_(0); }

        KOKKOS_INLINE_FUNCTION T& operator=(T other) { view_(0) = other; return view_(0);}

        KOKKOS_INLINE_FUNCTION T* data() { return view_.data(); }

        KOKKOS_INLINE_FUNCTION const T* data() const { return view_.data(); }

        KOKKOS_INLINE_FUNCTION View GetView() const { return *this; }

        KOKKOS_INLINE_FUNCTION operator ConstDeviceView() const { return { view_ }; }
    };

    KokkosScalarView(T init = T{})
        : data_(init),
          host_view_(&data_, 1)
    {
    }

    KOKKOS_INLINE_FUNCTION T* data()
    {
        return &data_;
    }

    KOKKOS_INLINE_FUNCTION const T* data() const
    {
        return &data_;
    }

    KOKKOS_INLINE_FUNCTION operator T() const { return data_; }

    constexpr KOKKOS_INLINE_FUNCTION Kokkos::View<T*> GetDeviceView() const { return device_view_; }

    void CopyToHost() const
    {
        Kokkos::deep_copy(host_view_, device_view_);
    }

    void CopyToDevice() const
    {
        Kokkos::deep_copy(device_view_, host_view_);
    }

    T& host_value() const { return data_; }
  };
}