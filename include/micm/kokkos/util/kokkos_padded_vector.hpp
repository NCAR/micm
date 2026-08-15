// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>
#include <micm/util/view_category.hpp>

#include <Kokkos_Core.hpp>
#include <initializer_list>
#include <vector>

namespace micm
{
  template<class T, Index L>
  class KokkosPaddedVector
  {
    mutable std::vector<T> data_;  // host data
    mutable Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> host_view_;
    mutable Kokkos::View<T*> device_view_;  // device data
    Index size_;

   public:
    using value_type = T;
    using category = PaddedVectorTag;
    struct DeviceView;
    struct ConstDeviceView;
    using ViewType = DeviceView;
    using ConstViewType = ConstDeviceView;

    struct ConstDeviceView
    {
      Kokkos::View<const T*> view_;
      Index size_;

      KOKKOS_INLINE_FUNCTION const T& operator[](Index i) const
      {
        return view_(i);
      }

      KOKKOS_INLINE_FUNCTION Index size() const  // NOLINT(readability-identifier-naming)
      {
        return size_;
      }

      // So the KokkosVectorLike concept accepts a DeviceView also
      KOKKOS_INLINE_FUNCTION const ConstDeviceView& GetView() const
      {
        return *this;
      }

      KOKKOS_INLINE_FUNCTION const T* begin() const  // NOLINT(readability-identifier-naming)
      {
        return view_.data();
      }
      KOKKOS_INLINE_FUNCTION const T* end() const  // NOLINT(readability-identifier-naming)
      {
        return view_.data() + size_;
      }
      KOKKOS_INLINE_FUNCTION const T* data() const  // NOLINT(readability-identifier-naming)
      {
        return view_.data();
      }
    };

    struct DeviceView
    {
      Kokkos::View<T*> view_;
      Index size_;

      KOKKOS_INLINE_FUNCTION T& operator[](Index i) const
      {
        return view_(i);
      }

      KOKKOS_INLINE_FUNCTION Index size() const  // NOLINT(readability-identifier-naming)
      {
        return size_;
      }

      // So the KokkosVectorLike concept accepts a DeviceView also
      KOKKOS_INLINE_FUNCTION DeviceView& GetView() const
      {
        return *this;
      }

      // Implicit conversion to ConstDeviceView (mirrors Kokkos::View<T*> -> View<const T*>)
      KOKKOS_INLINE_FUNCTION operator ConstDeviceView() const
      {
        return { view_ };
      }

      KOKKOS_INLINE_FUNCTION const T* begin() const  // NOLINT(readability-identifier-naming)
      {
        return view_.data();
      }
      KOKKOS_INLINE_FUNCTION const T* end() const  // NOLINT(readability-identifier-naming)
      {
        return view_.data() + size_;
      }
      KOKKOS_INLINE_FUNCTION const T* data() const  // NOLINT(readability-identifier-naming)
      {
        return view_.data();
      }
    };

    static constexpr Index GroupVectorSize()
    {
      return L;
    }

    KokkosPaddedVector()
        : size_(0)
    {
    }

    KokkosPaddedVector(Index n, T init = T{})
        : data_((((n + L - 1) / L) * L), init),
          host_view_(this->data_.data(), this->data_.size()),
          device_view_("padded_vector", this->data_.size()),
          size_(n)
    {
    }

    KokkosPaddedVector(std::initializer_list<T> init)
        : data_(((init.size() + L - 1) / L) * L, T{}),
          host_view_(this->data_.data(), this->data_.size()),
          device_view_("padded_vector", this->data_.size()),
          size_(init.size())
    {
      std::copy(init.begin(), init.end(), data_.begin());
    }

    KokkosPaddedVector(std::vector<T> init)
        : data_(((init.size() + L - 1) / L) * L, T{}),
          host_view_(this->data_.data(), this->data_.size()),
          device_view_("padded_vector", this->data_.size()),
          size_(init.size())
    {
      std::copy(init.begin(), init.end(), data_.begin());
    }

    Index size() const  // NOLINT(readability-identifier-naming)
    {
      return size_;
    }

    Index PaddedSize() const
    {
      return data_.size();
    }

    T& operator[](Index i)
    {
      return data_[i];
    }

    const T& operator[](Index i) const
    {
      return data_[i];
    }

    auto begin()  // NOLINT(readability-identifier-naming)
    {
      return data_.begin();
    }

    auto end()  // NOLINT(readability-identifier-naming)
    {
      return data_.begin() + size_;
    }

    auto begin() const  // NOLINT(readability-identifier-naming)
    {
      return data_.begin();
    }

    auto end() const  // NOLINT(readability-identifier-naming)
    {
      return data_.begin() + size_;
    }

    T* data()  // NOLINT(readability-identifier-naming)
    {
      return data_.data();
    }

    const T* data() const  // NOLINT(readability-identifier-naming)
    {
      return data_.data();
    }

    bool operator==(const KokkosPaddedVector<T, L>& other) const
    {
      return data_ == other.data_;
    }

    bool operator==(const std::vector<T>& other) const
    {
      return data_ == other;
    }

    friend bool operator==(const std::vector<T>& lhs, const KokkosPaddedVector<T, L>& rhs)
    {
      return rhs.data_ == lhs;
    }

    /// @brief Copy host data to the device view
    void CopyToDevice() const
    {
      Kokkos::deep_copy(device_view_, host_view_);
    }

    /// @brief Copy device data to the host vector
    void CopyToHost() const
    {
      Kokkos::deep_copy(host_view_, device_view_);
    }

    KOKKOS_INLINE_FUNCTION DeviceView GetView()
    {
      return { device_view_, size_ };
    }

    KOKKOS_INLINE_FUNCTION ConstDeviceView GetView() const
    {
      return { device_view_, size_ };
    }
  };
}  // namespace micm