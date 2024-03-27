// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <algorithm>
#include <cassert>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace micm
{
  /// Concept for vectorizable matrices
  template<typename T>
  concept VectorizableSparse = requires(T t)
  {
    t.GroupSize(0);
    t.GroupVectorSize();
    t.NumberOfGroups(0);
  };

  template<class T, class OrderingPolicy>
  class SparseMatrixBuilder;

  template<class T, class OrderingPolicy = SparseMatrixStandardOrdering>
  class SparseMatrix;

  template<class T>
  using StandardSparseMatrix = SparseMatrix<T, SparseMatrixStandardOrdering>;

  /// @brief A sparse block-diagonal 2D matrix class with contiguous memory
  ///
  /// Each block sub-matrix is square and has the same structure of non-zero elements
  ///
  /// Sparse matrix data structure follows the Compressed Sparse Row (CSR) pattern
  ///
  /// The template parameters are the type of the matrix elements and a class that
  /// defines the sizing and ordering of the data elements
  template<class T, class OrderingPolicy>
  class SparseMatrix : public OrderingPolicy
  {
   protected:
    std::size_t number_of_blocks_;        // Number of block sub-matrices in the overall matrix
    std::vector<std::size_t> row_ids_;    // Row indices of each non-zero element in a block
    std::vector<std::size_t> row_start_;  // Index in data_ and row_ids_ of the start of each row in a block
    std::vector<T> data_;                 // Value of each non-zero matrix element

   private:
    friend class SparseMatrixBuilder<T, OrderingPolicy>;
    friend class ProxyRow;
    friend class ConstProxyRow;
    friend class Proxy;
    friend class ConstProxy;

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

    class ConstProxy
    {
      const SparseMatrix& matrix_;
      std::size_t block_id_;
      std::size_t row_id_;

     public:
      ConstProxy(const SparseMatrix& matrix, std::size_t block_id, std::size_t row_id)
          : matrix_(matrix),
            block_id_(block_id),
            row_id_(row_id)
      {
      }

      std::size_t size() const
      {
        return matrix_.row_start_.size() - 1;
      }

      const T& operator[](std::size_t y) const
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

    class ConstProxyRow
    {
      const SparseMatrix& matrix_;
      std::size_t block_id_;

     public:
      ConstProxyRow(const SparseMatrix& matrix, std::size_t block_id)
          : matrix_(matrix),
            block_id_(block_id)
      {
      }

      std::size_t size() const
      {
        return matrix_.row_start_.size() - 1;
      }

      ConstProxy operator[](std::size_t x) const
      {
        return ConstProxy(matrix_, block_id_, x);
      }
    };

   public:
    static SparseMatrixBuilder<T, OrderingPolicy> create(std::size_t block_size)
    {
      return SparseMatrixBuilder<T, OrderingPolicy>{ block_size };
    }

    SparseMatrix() = default;

    SparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
        : number_of_blocks_(builder.number_of_blocks_),
          row_ids_(builder.RowIdsVector()),
          row_start_(builder.RowStartVector()),
          data_(OrderingPolicy::VectorSize(number_of_blocks_, row_ids_, row_start_), builder.initial_value_)
    {
    }

    SparseMatrix<T, OrderingPolicy>& operator=(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
    {
      number_of_blocks_ = builder.number_of_blocks_;
      row_ids_ = builder.RowIdsVector();
      row_start_ = builder.RowStartVector();
      data_ = std::vector<T>(OrderingPolicy::VectorSize(number_of_blocks_, row_ids_, row_start_), builder.initial_value_);

      return *this;
    }

    std::vector<T>& AsVector()
    {
      return data_;
    }

    const std::vector<T>& AsVector() const
    {
      return data_;
    }

    std::size_t VectorIndex(std::size_t block, std::size_t row, std::size_t column) const
    {
      return OrderingPolicy::VectorIndex(number_of_blocks_, row_ids_, row_start_, block, row, column);
    }

    std::size_t VectorIndex(std::size_t row, std::size_t column) const
    {
      if (number_of_blocks_ != 1)
        throw std::invalid_argument("Multi-block SparseMatrix access must specify block index");
      return VectorIndex(0, row, column);
    }

    bool IsZero(std::size_t row, std::size_t column) const
    {
      if (row >= row_start_.size() - 1 || column >= row_start_.size() - 1)
        throw std::invalid_argument("SparseMatrix element out of range");
      auto begin = std::next(row_ids_.begin(), row_start_[row]);
      auto end = std::next(row_ids_.begin(), row_start_[row + 1]);
      auto elem = std::find(begin, end, column);
      if (elem == end)
        return true;
      return false;
    }

    std::size_t size() const
    {
      return number_of_blocks_;
    }

    std::size_t FlatBlockSize() const
    {
      return row_ids_.size();
    }

    ConstProxyRow operator[](std::size_t b) const
    {
      return ConstProxyRow(*this, b);
    }

    ProxyRow operator[](std::size_t b)
    {
      return ProxyRow(*this, b);
    }

    SparseMatrix& operator=(T val)
    {
      std::transform(data_.begin(), data_.end(), data_.begin(), [&](auto& _) { return val; });
      return *this;
    }

    const std::vector<std::size_t>& RowStartVector() const
    {
      return row_start_;
    }

    const std::vector<std::size_t>& RowIdsVector() const
    {
      return row_ids_;
    }
  };

  template<class T, class OrderingPolicy = SparseMatrixStandardOrdering>
  class SparseMatrixBuilder
  {
    std::size_t number_of_blocks_{ 1 };
    std::size_t block_size_;
    std::set<std::pair<std::size_t, std::size_t>> non_zero_elements_{};
    T initial_value_{};
    friend class SparseMatrix<T, OrderingPolicy>;

   public:
    SparseMatrixBuilder() = delete;

    SparseMatrixBuilder(std::size_t block_size)
        : block_size_(block_size)
    {
    }

    operator SparseMatrix<T, OrderingPolicy>() const
    {
      return SparseMatrix<T, OrderingPolicy>(*this);
    }

    SparseMatrixBuilder& number_of_blocks(std::size_t n)
    {
      number_of_blocks_ = n;
      return *this;
    }

    SparseMatrixBuilder& with_element(std::size_t x, std::size_t y)
    {
      if (x >= block_size_ || y >= block_size_)
        throw std::invalid_argument("SparseMatrix element out of range");
      non_zero_elements_.insert(std::make_pair(x, y));
      return *this;
    }

    SparseMatrixBuilder& initial_value(T inital_value)
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
