// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cassert>
#include <cmath>
#include <vector>

#include <micm/util/exit_codes.hpp>

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
        if (other.size() < y_dim_) {
          std::cerr << "Matrix row size mismatch in assignment from vector";
          std::exit(micm::ExitCodes::InvalidMatrixDimension);
        }
        auto iter = std::next(matrix_.data_.begin(), block_index_ * y_dim_ * L + row_index_);
#if __cpluscplus >= 201907L
	std::for_each(other.begin(), std::next(other.begin(), y_dim_), [&](T const& elem) {
          *iter = elem;
          iter += L;
        });
#else
	for(auto elem = other.begin(); elem < std::next(other.begin(), y_dim_); ++elem)
	{
	  *iter = *elem;
	  iter += L;
	}
#endif
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
                  if (other_row.size() != y_dim)
                  {
                    std::cerr << "Invalid vector for matrix assignment\n";
                    std::exit(micm::ExitCodes::InvalidMatrixDimension);
                  }
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

    std::size_t NumberOfBlocks() const
    {
      return std::ceil(x_dim_ / (double)L);
    }

    std::size_t BlockSize() const
    {
      return L * y_dim_;
    }

    constexpr std::size_t VectorSize() const
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
