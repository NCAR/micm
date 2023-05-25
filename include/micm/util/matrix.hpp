// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <vector>

namespace micm
{

  /// @brief A 2D array class with contiguous memory
  template<class T>
  class Matrix
  {
    std::vector<T> data_;
    const std::size_t x_dim_;
    const std::size_t y_dim_;

    friend class Proxy;

    class Proxy
    {
      Matrix &matrix_;
      std::size_t offset_;

     public:
      Proxy(Matrix &matrix, std::size_t offset)
          : matrix_(matrix),
            offset_(offset)
      {
      }
      T &operator[](std::size_t y)
      {
        return matrix_.data_[offset_ + y];
      }
    };

   public:
    Matrix(std::size_t x_dim, std::size_t y_dim)
        : x_dim_(x_dim),
          y_dim_(y_dim),
          data_(x_dim * y_dim)
    {
    }
    Proxy operator[](std::size_t x)
    {
      return Proxy(*this, x * y_dim_);
    }

    std::vector<T> &AsVector()
    {
      return data_;
    }
  };

}  // namespace micm