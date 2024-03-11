// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#ifndef MICM_DEFAULT_VECTOR_SIZE
#  define MICM_DEFAULT_VECTOR_SIZE 4
#endif

namespace micm
{

  /// @brief A 2D array class with contiguous memory structured to encourage vectorization
  ///
  /// The memory layout groups rows into groups whose size can be set such that for a single
  /// column, the group of rows can fit in the vector register.
  ///
  /// The template arguments are the type of the matrix elements and the size of the number
  /// of rows per group.
  template<class T, std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  class VectorMatrix
  {
   protected:
    std::vector<T> data_;
    std::size_t x_dim_;  // number of rows
    std::size_t y_dim_;  // number of columns

   private:
    friend class Proxy;
    friend class ConstProxy;

    class Proxy
    {
      VectorMatrix &matrix_;
      std::size_t group_index_;
      std::size_t row_index_;
      std::size_t y_dim_;

     public:
      Proxy(VectorMatrix &matrix, std::size_t group_index, std::size_t row_index, std::size_t y_dim)
          : matrix_(matrix),
            group_index_(group_index),
            row_index_(row_index),
            y_dim_(y_dim)
      {
      }
      Proxy &operator=(const std::vector<T> other)
      {
        if (other.size() < y_dim_)
        {
          throw std::runtime_error("Matrix row size mismatch in assignment from vector.");
        }
        auto iter = std::next(matrix_.data_.begin(), group_index_ * y_dim_ * L + row_index_);
        std::for_each(
            other.begin(),
            std::next(other.begin(), y_dim_),
            [&](T const &elem)
            {
              *iter = elem;
              // don't iterate past the end of the vector
              std::size_t remaining_elements = std::distance(iter, matrix_.data_.end());
              iter += std::min(L, remaining_elements);
            });
        return *this;
      }
      operator std::vector<T>() const
      {
        std::vector<T> vec(y_dim_);
        auto iter = std::next(matrix_.data_.begin(), group_index_ * y_dim_ * L + row_index_);
        for (auto &elem : vec)
        {
          elem = *iter;
          // don't iterate past the end of the vector
          std::size_t remaining_elements = std::distance(iter, matrix_.data_.end());
          iter += std::min(L, remaining_elements);
        }
        return vec;
      }
      std::size_t size() const
      {
        return y_dim_;
      }
      T &operator[](std::size_t y)
      {
        return matrix_.data_[(group_index_ * y_dim_ + y) * L + row_index_];
      }
    };

    class ConstProxy
    {
      const VectorMatrix &matrix_;
      std::size_t group_index_;
      std::size_t row_index_;
      std::size_t y_dim_;

     public:
      ConstProxy(const VectorMatrix &matrix, std::size_t group_index, std::size_t row_index, std::size_t y_dim)
          : matrix_(matrix),
            group_index_(group_index),
            row_index_(row_index),
            y_dim_(y_dim)
      {
      }
      operator std::vector<T>() const
      {
        std::vector<T> vec(y_dim_);
        auto iter = std::next(matrix_.data_.begin(), group_index_ * y_dim_ * L + row_index_);
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
        return matrix_.data_[(group_index_ * y_dim_ + y) * L + row_index_];
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
                  if (other_row.size() != y_dim)
                  {
                    throw std::runtime_error("Invalid vector for matrix assignment");
                  }
                  auto iter = std::next(data.begin(), std::floor(i_row / (double)L) * y_dim * L + i_row % L);
                  for (auto &elem : other_row)
                  {
                    *iter = elem;
                    // don't iterate past the end of the vector
                    std::size_t remaining_elements = std::distance(iter, data.end());
                    iter += std::min(L, remaining_elements);
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

    std::size_t NumberOfGroups() const
    {
      return std::ceil(x_dim_ / (double)L);
    }

    std::size_t GroupSize() const
    {
      return L * y_dim_;
    }

    constexpr std::size_t GroupVectorSize() const
    {
      return L;
    }

    ConstProxy operator[](std::size_t x) const
    {
      return ConstProxy(*this, std::floor(x / L), x % L, y_dim_);
    }

    Proxy operator[](std::size_t x)
    {
      return Proxy(*this, std::floor(x / L), x % L, y_dim_);
    }

    VectorMatrix &operator=(T val)
    {
      std::transform(data_.begin(), data_.end(), data_.begin(), [&](auto &_) { return val; });
      return *this;
    }

    void ForEach(const std::function<void(T &, const T &)> f, const VectorMatrix &a)
    {
      auto this_iter = data_.begin();
      auto a_iter = a.AsVector().begin();
      const std::size_t n = std::floor(x_dim_ / L) * L * y_dim_;
      for (std::size_t i = 0; i < n; ++i)
        f(*(this_iter++), *(a_iter++));
      const std::size_t l = x_dim_ % L;
      for (std::size_t y = 0; y < y_dim_; ++y)
        for (std::size_t x = 0; x < l; ++x)
          f(this_iter[y * L + x], a_iter[y * L + x]);
    }

    void ForEach(const std::function<void(T &, const T &, const T &)> f, const VectorMatrix &a, const VectorMatrix &b)
    {
      auto this_iter = data_.begin();
      auto a_iter = a.AsVector().begin();
      auto b_iter = b.AsVector().begin();
      const std::size_t n = std::floor(x_dim_ / L) * L * y_dim_;
      for (std::size_t i = 0; i < n; ++i)
        f(*(this_iter++), *(a_iter++), *(b_iter++));
      const std::size_t l = x_dim_ % L;
      for (std::size_t y = 0; y < y_dim_; ++y)
        for (std::size_t x = 0; x < l; ++x)
          f(this_iter[y * L + x], a_iter[y * L + x], b_iter[y * L + x]);
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
