// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <vector>

namespace micm
{

  /// @brief A 2D array class with contiguous memory
  template<class T, std::size_t x_dim, std::size_t y_dim>
  class Matrix
  {
    std::vector<T> data_;

    friend class Proxy;

    class Proxy
    {
      Matrix &matrix_;
      std::size_t x_;

     public:
      Proxy(Matrix &matrix, std::size_t x)
          : matrix_(matrix),
            x_(x)
      {
      }
      T &operator[](std::size_t y)
      {
        return matrix_.data_[x_ * y_dim + y];
      }
    };

   public:
    Matrix()
        : data_(x_dim * y_dim)
    {
    }
    Proxy operator[](std::size_t x)
    {
      return Proxy(*this, x);
    }

    std::vector<T> &AsVector()
    {
      return data_;
    }
  };

}  // namespace micm