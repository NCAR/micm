// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cassert>
#include <stdexcept>
#include <utility>
#include <vector>

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

   public:
    static SparseMatrixBuilder<T> create(std::size_t block_size)
    {
      return SparseMatrixBuilder<T>{ block_size };
    }

    SparseMatrix(SparseMatrixBuilder<T>& builder)
        : number_of_blocks_(builder.number_of_blocks_),
          data_(builder.NumberOfElements(), builder.initial_value_),
          row_ids_(builder.RowIdsVector()),
          row_start_(builder.RowStartVector())
    {
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