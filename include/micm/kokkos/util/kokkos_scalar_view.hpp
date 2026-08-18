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
    mutable Kokkos::View<T*, Kokkos::HostSpace> data_{ "scalar_host", 1 };
    mutable Kokkos::View<T*> device_view_{ "scalar", 1 };

   public:
    using value_type = T;
    template<class U>
    struct View;
    using ViewType = View<T>;
    using ConstViewType = View<const T>;

    KokkosScalarView(T init = T{})
    {
      data_(0) = init;
    }

    KOKKOS_INLINE_FUNCTION T* data()  // NOLINT(readability-identifier-naming)
    {
      KOKKOS_IF_ON_DEVICE((return device_view_.data();))
      KOKKOS_IF_ON_HOST((return data_.data();))
    }

    KOKKOS_INLINE_FUNCTION const T* data() const  // NOLINT(readability-identifier-naming)
    {
      KOKKOS_IF_ON_DEVICE((return device_view_.data();))
      KOKKOS_IF_ON_HOST((return data_.data();))
    }

    KOKKOS_INLINE_FUNCTION operator T() const
    {
      KOKKOS_IF_ON_DEVICE((return device_view_(0);))
      KOKKOS_IF_ON_HOST((return data_(0);))
    }

    T& operator=(const T& other)
    {
      data_(0) = other;
      return data_(0);
    }

    KokkosScalarView<T>& operator=(const KokkosScalarView<T>& other)
    {
      if (this != &other)
      {
        data_(0) = other.data_(0);
      }
      return *this;
    }

    constexpr KOKKOS_INLINE_FUNCTION Kokkos::View<T*> GetDeviceView() const
    {
      return device_view_;
    }

    void CopyToHost() const
    {
      Kokkos::deep_copy(data_, device_view_);
    }

    void CopyToDevice() const
    {
      Kokkos::deep_copy(device_view_, data_);
    }

    T& HostValue() const
    {
      return data_(0);
    }
  };
}  // namespace micm