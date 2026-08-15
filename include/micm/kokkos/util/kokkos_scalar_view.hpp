// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

// NOLINTNEXTLINE(clang-diagnostic-error): Kokkos isn't included in the clang-tidy build
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
    struct View;
    struct ConstView;
    using ViewType = View;
    using ConstViewType = ConstView;

    struct ConstDeviceView
    {
      Kokkos::View<const T*> view_;

      KOKKOS_INLINE_FUNCTION operator T() const
      {
        return view_(0);
      }

      KOKKOS_INLINE_FUNCTION const T* data() const  // NOLINT(readability-identifier-naming)
      {
        return view_.data();
      }

      KOKKOS_INLINE_FUNCTION ConstView GetView() const
      {
        return *this;
      }
    };

    struct DeviceView
    {
      Kokkos::View<T*> view_;

      KOKKOS_INLINE_FUNCTION operator T() const
      {
        return view_(0);
      }

      KOKKOS_INLINE_FUNCTION T& operator=(const T& other)
      {
        view_(0) = other;
        return view_(0);
      }

      KOKKOS_INLINE_FUNCTION DeviceView& operator=(const DeviceView& other)
      {
        view_(0) = other.view_(0);
        return *this;
      }

      KOKKOS_INLINE_FUNCTION T* data()  // NOLINT(readability-identifier-naming)
      {
        return view_.data();
      }

      KOKKOS_INLINE_FUNCTION const T* data() const  // NOLINT(readability-identifier-naming)
      {
        return view_.data();
      }

      KOKKOS_INLINE_FUNCTION View GetView() const
      {
        return *this;
      }

      KOKKOS_INLINE_FUNCTION operator ConstDeviceView() const
      {
        return { view_ };
      }
    };

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
      data_(0) = other.data_(0);
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