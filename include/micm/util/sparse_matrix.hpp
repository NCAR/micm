// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <algorithm>
#include <cassert>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include <set>

namespace micm
{

  template<class T>
  class SparseMatrixBuilder;

  /// @brief A sparse block-diagonal 2D matrix class with contiguous memory
  ///
  /// Each block sub-matrix is square and has the same structure of non-zero elements
  ///
  /// Sparse matrix data structure follows the Compressed Sparse Row (CSR) pattern
  template<class T>
  class SparseMatrix
  {
    std::size_t number_of_blocks_;        // Number of block sub-matrices in the overall matrix
    std::vector<T> data_;                 // Value of each non-zero matrix element
    std::vector<std::size_t> row_ids_;    // Row indices of each non-zero element in a block
    std::vector<std::size_t> row_start_;  // Index in data_ and row_ids_ of the start of each column in a block

    friend class SparseMatrixBuilder<T>;
    friend class ProxyRow;
    friend class Proxy;

    class Proxy
    {
      SparseMatrix& matrix_;
      std::size_t block_id_;
      std::size_t row_id_;

     public:
      Proxy(SparseMatrix& matrix, std::size_t block_id, std::size_t row_id)
          : matrix_(matrix),
            block_id_(block_id),
            row_id_(row_id)
      {
      }

      std::size_t size() const
      {
        return matrix_.row_start_.size() - 1;
      }

      T& operator[](std::size_t y)
      {
        return matrix_.data_[matrix_.VectorIndex(block_id_, row_id_, y)];
      }
    };

    class ProxyRow
    {
      SparseMatrix& matrix_;
      std::size_t block_id_;

     public:
      ProxyRow(SparseMatrix& matrix, std::size_t block_id)
          : matrix_(matrix),
            block_id_(block_id)
      {
      }

      std::size_t size() const
      {
        return matrix_.row_start_.size() - 1;
      }

      Proxy operator[](std::size_t x)
      {
        return Proxy(matrix_, block_id_, x);
      }
    };

   public:
    static SparseMatrixBuilder<T> create(std::size_t block_size)
    {
      return SparseMatrixBuilder<T>{ block_size };
    }

    SparseMatrix() = default;

    SparseMatrix(SparseMatrixBuilder<T>& builder)
        : number_of_blocks_(builder.number_of_blocks_),
          data_(builder.NumberOfElements(), builder.initial_value_),
          row_ids_(builder.RowIdsVector()),
          row_start_(builder.RowStartVector())
    {
    }

    std::vector<T>& AsVector()
    {
      return data_;
    }

    std::size_t VectorIndex(std::size_t block, std::size_t row, std::size_t column) const
    {
      if (row >= row_start_.size() - 1 || column >= row_start_.size() - 1 || block >= number_of_blocks_)
        throw std::invalid_argument("SparseMatrix element out of range");
      auto begin = std::next(row_ids_.begin(), row_start_[row]);
      auto end = std::next(row_ids_.begin(), row_start_[row + 1]);
      auto elem = std::find(begin, end, column);
      if (elem == end)
        throw std::invalid_argument("SparseMatrix zero element access not allowed");
      return std::size_t{ (elem - row_ids_.begin()) + block * row_ids_.size() };
    }

    std::size_t VectorIndex(std::size_t row, std::size_t column) const
    {
      if (number_of_blocks_ != 1)
        throw std::invalid_argument("Multi-block SparseMatrix access must specify block index");
      return VectorIndex(0, row, column);
    }

    std::size_t size() const
    {
      return number_of_blocks_;
    }

    std::size_t FlatBlockSize() const
    {
      return row_ids_.size();
    }

    ProxyRow operator[](std::size_t b)
    {
      return ProxyRow(*this, b);
    }
  };

  template<class T>
  class SparseMatrixBuilder
  {
    std::size_t number_of_blocks_{ 1 };
    std::size_t block_size_;
    std::set<std::pair<std::size_t, std::size_t>> non_zero_elements_{};
    T initial_value_{};
    friend class SparseMatrix<T>;

   public:
    SparseMatrixBuilder() = delete;

    SparseMatrixBuilder(std::size_t block_size)
        : block_size_(block_size)
    {
    }

    operator SparseMatrix<T>() const
    {
      return SparseMatrix<T>(*this);
    }

    SparseMatrixBuilder<T>& number_of_blocks(std::size_t n)
    {
      number_of_blocks_ = n;
      return *this;
    }

    SparseMatrixBuilder<T>& with_element(std::size_t x, std::size_t y)
    {
      if (x >= block_size_ || y >= block_size_)
        throw std::invalid_argument("SparseMatrix element out of range");
      non_zero_elements_.insert(std::make_pair(x, y));
      return *this;
    }

    SparseMatrixBuilder<T>& initial_value(T inital_value)
    {
      initial_value_ = inital_value;
      return *this;
    }

    std::size_t NumberOfElements() const
    {
      return non_zero_elements_.size() * number_of_blocks_;
    }

    std::vector<std::size_t> RowIdsVector() const
    {
      std::vector<std::size_t> ids;
      ids.reserve(non_zero_elements_.size());
      std::transform(
          non_zero_elements_.begin(),
          non_zero_elements_.end(),
          std::back_inserter(ids),
          [](const std::pair<std::size_t, std::size_t>& elem) { return elem.second; });
      return ids;
    }
    std::vector<std::size_t> RowStartVector() const
    {
      std::vector<std::size_t> starts(block_size_ + 1, 0);
      std::size_t total_elem = 0;
      std::size_t curr_row = 0;
      for (auto& elem : non_zero_elements_)
      {
        while (curr_row < elem.first)
          starts[(curr_row++) + 1] = total_elem;
        ++total_elem;
      }
      starts[curr_row + 1] = total_elem;
      return starts;
    }
  };

}  // namespace micm