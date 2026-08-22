// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>
#include <micm/util/view_category.hpp>

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <initializer_list>
#include <type_traits>
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
    template<class U>
    struct DeviceView;
    using ViewType = DeviceView<T>;
    using ConstViewType = DeviceView<const T>;

    template<class U>
    struct DeviceView
    {
      Kokkos::View<U*> view_;
      Index size_;

      KOKKOS_INLINE_FUNCTION U& operator[](Index i) const  // NOLINT(readability-identifier-naming)
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

      // NOLINTNEXTLINE(modernize-use-constraints) nvhpc warnings when constraints are used
      template<class V = U, std::enable_if_t<!std::is_const_v<V>, int> = 0>
      KOKKOS_INLINE_FUNCTION operator DeviceView<const U>() const
      {
        // size_ has to be carried across explicitly. DeviceView is an aggregate, so a
        // brace-init that names only view_ leaves size_ zero-initialized and the converted
        // view reports size() == 0 while still pointing at the full padded allocation.
        return { view_, size_ };
      }

      KOKKOS_INLINE_FUNCTION const U* begin() const  // NOLINT(readability-identifier-naming)
      {
        return view_.data();
      }
      KOKKOS_INLINE_FUNCTION const U* end() const  // NOLINT(readability-identifier-naming)
      {
        return view_.data() + size_;
      }
      KOKKOS_INLINE_FUNCTION const U* data() const  // NOLINT(readability-identifier-naming)
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

    KokkosPaddedVector(const KokkosPaddedVector& other)
        : data_(other.data_),
          host_view_(this->data_.data(), this->data_.size()),
          device_view_("padded_vector", this->data_.size()),
          size_(other.size_)
    {
      Kokkos::deep_copy(device_view_, other.device_view_);
    }

    KokkosPaddedVector& operator=(const KokkosPaddedVector& other)
    {
      if (this == &other)
      {
        return *this;
      }
      data_ = other.data_;
      // size_ is the logical length and is independent of data_.size(), which is padded.
      // Leaving it behind here would silently keep the assignee's old length -- State's copy
      // assignment reaches this operator through conditions_ and the diagonal-element
      // vectors, so a stale size_ corrupts every kernel that bounds its loop by size().
      size_ = other.size_;
      host_view_ = Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
          this->data_.data(), this->data_.size());
      Kokkos::realloc(device_view_, other.device_view_.extent(0));
      Kokkos::deep_copy(device_view_, other.device_view_);
      return *this;
    }

    KokkosPaddedVector(KokkosPaddedVector&& other) = default;

    KokkosPaddedVector& operator=(KokkosPaddedVector&& other) = default;

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

    KOKKOS_INLINE_FUNCTION DeviceView<T> GetView()
    {
      return { device_view_, size_ };
    }

    KOKKOS_INLINE_FUNCTION DeviceView<const T> GetView() const
    {
      return { device_view_, size_ };
    }
  };
}  // namespace micm