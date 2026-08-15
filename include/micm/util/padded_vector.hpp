// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>
#include <micm/util/view_category.hpp>

#include <initializer_list>
#include <vector>

namespace micm
{
  /// @brief A vector class with padded cells for use in Matrix::Function lambdas
  template<class T, Index L>
  class PaddedVector
  {
    std::vector<T> data_;
    Index size_;

   public:
    using value_type = T;
    struct View;
    struct ConstView;
    using ViewType = View;
    using ConstViewType = ConstView;
    using category = PaddedVectorTag;

    struct ConstView
    {
      const T* data_;
      Index size_;
      Index padded_size_;

      const T& operator[](Index i) const
      {
        return data_[i];
      }
      Index size() const  // NOLINT(readability-identifier-naming)
      {
        return size_;
      }
      Index PaddedSize() const
      {
        return padded_size_;
      }
      ConstView GetView() const
      {
        return *this;
      }

      const T* begin() const  // NOLINT(readability-identifier-naming)
      {
        return data_;
      }
      const T* end() const  // NOLINT(readability-identifier-naming)
      {
        return data_ + size_;
      }
      const T* data() const  // NOLINT(readability-identifier-naming)
      {
        return data_;
      }
    };

    struct View
    {
      T* data_;
      Index size_;
      Index padded_size_;

      T& operator[](Index i)
      {
        return data_[i];
      }
      const T& operator[](Index i) const
      {
        return data_[i];
      }
      Index size() const  // NOLINT(readability-identifier-naming)
      {
        return size_;
      }
      Index PaddedSize() const
      {
        return padded_size_;
      }
      View GetView() const
      {
        return *this;
      }

      T* begin()  // NOLINT(readability-identifier-naming)
      {
        return data_;
      }
      T* end()  // NOLINT(readability-identifier-naming)
      {
        return data_ + size_;
      }
      const T* begin() const  // NOLINT(readability-identifier-naming)
      {
        return data_;
      }
      const T* end() const  // NOLINT(readability-identifier-naming)
      {
        return data_ + size_;
      }
      T* data()  // NOLINT(readability-identifier-naming)
      {
        return data_;
      }
      const T* data() const  // NOLINT(readability-identifier-naming)
      {
        return data_;
      }

      // Implicit conversion to ConstView (mirrors DeviceView -> ConstDeviceView)
      operator ConstView() const
      {
        return { data_, size_, padded_size_ };
      }
    };

    static constexpr Index GroupVectorSize()
    {
      return L;
    }

    PaddedVector()
        : size_(0)
    {
    }

    PaddedVector(Index n, T init = T{})
        : data_(((n + L - 1) / L) * L, init),
          size_(n)
    {
    }

    PaddedVector(std::initializer_list<T> init)
        : data_(((init.size() + L - 1) / L) * L, T{}),
          size_(init.size())
    {
      std::copy(init.begin(), init.end(), data_.begin());
    }

    PaddedVector(const std::vector<T>& init)
        : data_(((init.size() + L - 1) / L) * L, T{}),
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

    bool empty() const  // NOLINT(readability-identifier-naming)
    {
      return data_.empty();
    }

    bool operator==(const PaddedVector<T, L>& other) const
    {
      return data_ == other.data_;
    }

    bool operator==(const std::vector<T>& other) const
    {
      return data_ == other;
    }

    friend bool operator==(const std::vector<T>& lhs, const PaddedVector<T, L>& rhs)
    {
      return rhs.data_ == lhs;
    }

    View GetView()
    {
      return { data_.data(), size_, (Index)data_.size() };
    }
    ConstView GetView() const
    {
      return { data_.data(), size_, (Index)data_.size() };
    }

    void CopyToHost() const
    {
    }

    void CopyToDevice() const
    {
    }
  };
}  // namespace micm