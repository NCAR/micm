// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>
#include <micm/util/view_category.hpp>

#include <initializer_list>
#include <vector>

namespace micm
{
  /// @brief A vector class with padded cells for use in VectorMatrix::Function lambdas
  template<class T, Index L>
  class PaddedVector
  {
    std::vector<T> data_;
    Index size_;

   public:
    using value_type = T;
    using category = PaddedVectorTag;
    
    static constexpr Index GroupVectorSize()
    {
        return L;
    }

    PaddedVector()
        : size_(0)
    {
    }

    PaddedVector(Index n, T init = T{})
        : data_(((n + L - 1) / L) * L, init), size_(n) { }

    PaddedVector(std::initializer_list<T> init)
        : data_(((init.size() + L - 1) / L) * L, T{}),
          size_(init.size())
    {
        std::copy(init.begin(), init.end(), data_.begin());
    }

    Index size() const
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

    auto begin()
    {
        return data_.begin();
    }

    auto end()
    {
        return data_.begin() + size_;
    }

    auto begin() const
    {
        return data_.begin();
    }

    auto end() const
    {
        return data_.begin() + size_;
    }

    T* data()
    {
        return data_.data();
    }

    const T* data() const
    {
        return data_.data();
    }

    void CopyToHost() const
    {
    }

    void CopyToDevice() const
    {
    }
  };
}