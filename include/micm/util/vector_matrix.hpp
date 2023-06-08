// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cassert>
#include <cmath>
#include <vector>

namespace micm
{

  /// @brief A 2D array class with contiguous memory structured to encourage vectorization
  ///
  /// The memory layout groups rows into blocks whose size can be set such that for a single
  /// column, the block of rows can fit in the vector register.
  ///
  /// The template arguments are the type of the matrix elements and the size of the number
  /// of rows per block.
  template<class T, std::size_t L>
  class VectorMatrix
  {
    std::vector<T> data_;
    std::size_t x_dim_;
    std::size_t y_dim_;

    friend class Proxy;
    friend class ConstProxy;

    class Proxy
    {
      VectorMatrix &matrix_;
      std::size_t block_index_;
      std::size_t row_index_;
      std::size_t y_dim_;

     public:
      Proxy(VectorMatrix &matrix, std::size_t block_index, std::size_t row_index, std::size_t y_dim)
          : matrix_(matrix),
            block_index_(block_index),
            row_index_(row_index),
            y_dim_(y_dim)
      {
      }
      Proxy &operator=(const std::vector<T> other)
      {
        assert(other.size() >= y_dim_ && "VectorMatrix row size mismatch in assignment from vector");
        auto iter = std::next(matrix_.data_.begin(), block_index_ * y_dim_ * L + row_index_);
        for (auto &elem : other)
        {
          *iter = elem;
          iter += L;
        }
        return *this;
      }
      operator std::vector<T>() const
      {
        std::vector<T> vec(y_dim_);
        auto iter = std::next(matrix_.data_.begin(), block_index_ * y_dim_ * L + row_index_);
        for (auto &elem : vec)
        {
          elem = *iter;
          iter += L;
        }
        return vec;
      }
      std::size_t size() const
      {
        return y_dim_;
      }
      T &operator[](std::size_t y)
      {
        return matrix_.data_[(block_index_ * y_dim_ + y) * L + row_index_];
      }
    };

    class ConstProxy
    {
      const VectorMatrix &matrix_;
      std::size_t block_index_;
      std::size_t row_index_;
      std::size_t y_dim_;

     public:
      ConstProxy(const VectorMatrix &matrix, std::size_t block_index, std::size_t row_index, std::size_t y_dim)
          : matrix_(matrix),
            block_index_(block_index),
            row_index_(row_index),
            y_dim_(y_dim)
      {
      }
      operator std::vector<T>() const
      {
        std::vector<T> vec(y_dim_);
        auto iter = std::next(matrix_.data_.begin(), block_index_ * y_dim_ * L + row_index_);
        for (auto &elem : vec)
        {
          elem = *iter;
          iter += L;
        }
        return vec;
      }
      std::size_t size() const
      {
        return y_dim_;
      }
      const T &operator[](std::size_t y) const
      {
        return matrix_.data_[(block_index_ * y_dim_ + y) * L + row_index_];
      }
    };

   public:
    VectorMatrix()
        : x_dim_(0),
          y_dim_(0),
          data_()
    {
    }

    VectorMatrix(std::size_t x_dim, std::size_t y_dim)
        : x_dim_(x_dim),
          y_dim_(y_dim),
          data_(std::ceil(x_dim / (double)L) * L * y_dim)
    {
    }

    VectorMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value)
        : x_dim_(x_dim),
          y_dim_(y_dim),
          data_(std::ceil(x_dim / (double)L) * L * y_dim, initial_value)
    {
    }

    VectorMatrix(const std::vector<std::vector<T>> other)
        : x_dim_(other.size()),
          y_dim_(other.size() == 0 ? 0 : other[0].size()),
          data_(
              [&]() -> std::vector<T>
              {
                std::size_t x_dim = other.size();
                if (x_dim == 0)
                  return std::vector<T>(0);
                std::size_t y_dim = other[0].size();
                std::vector<T> data(std::ceil(x_dim / (double)L) * L * y_dim);
                std::size_t i_row = 0;
                for (auto &other_row : other)
                {
                  assert(other_row.size() == y_dim && "Invalid vector for matrix assignment");
                  auto iter = std::next(data.begin(), std::floor(i_row / (double)L) * y_dim * L + i_row % L);
                  for (auto &elem : other_row)
                  {
                    *iter = elem;
                    iter += L;
                  }
                  ++i_row;
                }
                return data;
              }())
    {
    }
    std::size_t size() const
    {
      return x_dim_;
    }

    ConstProxy operator[](std::size_t x) const
    {
      return ConstProxy(*this, std::floor(x / L), x % L, y_dim_);
    }

    Proxy operator[](std::size_t x)
    {
      return Proxy(*this, std::floor(x / L), x % L, y_dim_);
    }

    std::vector<T> &AsVector()
    {
      return data_;
    }

    const std::vector<T> &AsVector() const
    {
      return data_;
    }
  };

}  // namespace micm