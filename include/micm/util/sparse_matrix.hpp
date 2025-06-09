// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/matrix_error.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <ostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace micm
{
  /// Concept for vectorizable matrices
  template<typename T>
  concept VectorizableSparse = requires(T t) {
    t.GroupSize();
    t.GroupVectorSize();
    t.NumberOfGroups(0);
  };

  template<typename T>
  concept SparseMatrixConcept = requires(T t) {
    t.NumRows();
    t.NumColumns();
    t.NumberOfBlocks();
  };

  template<class T, class OrderingPolicy>
  class SparseMatrixBuilder;

  template<class T, class OrderingPolicy = SparseMatrixStandardOrdering>
  class SparseMatrix;

  using StandardSparseMatrix = SparseMatrix<double, SparseMatrixStandardOrdering>;

  /// @brief A sparse block-diagonal 2D matrix class with contiguous memory
  ///
  /// Each block sub-matrix is square and has the same structure of non-zero elements
  ///
  /// The template parameters are the type of the matrix elements and a class that
  /// defines the sizing and ordering of the data elements
  template<class T = double, class OrderingPolicy>
  class SparseMatrix : public OrderingPolicy
  {
   public:
    // Diagonal markowitz reordering requires an int argument, make sure one is always accessible
    using IntMatrix = SparseMatrix<int, OrderingPolicy>;
    using value_type = T;

   protected:
    std::size_t number_of_blocks_;  // Number of block sub-matrices in the overall matrix
    std::size_t block_size_;        // Size of each block sub-matrix (number of rows or columns per block)
    std::size_t number_of_non_zero_elements_per_block_;  // Number of non-zero elements in each block
    std::vector<T> data_;                                // Value of each non-zero matrix element

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

      std::size_t Size() const
      {
        return matrix_.block_size_;
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

      std::size_t Size() const
      {
        return matrix_.block_size_;
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

      std::size_t Size() const
      {
        return matrix_.block_size_;
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

      std::size_t Size() const
      {
        return matrix_.block_size_;
      }

      ConstProxy operator[](std::size_t x) const
      {
        return ConstProxy(matrix_, block_id_, x);
      }
    };

   public:
    static SparseMatrixBuilder<T, OrderingPolicy> Create(std::size_t block_size)
    {
      return SparseMatrixBuilder<T, OrderingPolicy>{ block_size };
    }

    SparseMatrix() = default;

    SparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
        : OrderingPolicy(builder.number_of_blocks_, builder.block_size_, builder.non_zero_elements_),
          number_of_blocks_(builder.number_of_blocks_),
          block_size_(builder.block_size_),
          number_of_non_zero_elements_per_block_(builder.non_zero_elements_.size()),
          data_(OrderingPolicy::VectorSize(number_of_blocks_), builder.initial_value_)
    {
    }

    SparseMatrix<T, OrderingPolicy>& operator=(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
    {
      OrderingPolicy::operator=(std::make_tuple(builder.number_of_blocks_, builder.block_size_, builder.non_zero_elements_));
      number_of_blocks_ = builder.number_of_blocks_;
      block_size_ = builder.block_size_;
      number_of_non_zero_elements_per_block_ = builder.non_zero_elements_.size();
      data_ = std::vector<T>(OrderingPolicy::VectorSize(number_of_blocks_), builder.initial_value_);

      return *this;
    }

    std::vector<std::size_t> DiagonalIndices(const std::size_t block_id) const
    {
      return OrderingPolicy::DiagonalIndices(number_of_blocks_, block_id);
    }

    void AddToDiagonal(T value)
    {
      OrderingPolicy::AddToDiagonal(number_of_blocks_, data_, value);
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
      return OrderingPolicy::VectorIndex(number_of_blocks_, block, row, column);
    }

    std::size_t VectorIndex(std::size_t row, std::size_t column) const
    {
      if (number_of_blocks_ != 1)
        throw std::system_error(make_error_code(MicmMatrixErrc::MissingBlockIndex));
      return VectorIndex(0, row, column);
    }

    std::size_t NumberOfBlocks() const
    {
      return number_of_blocks_;
    }

    std::size_t NumRows() const
    {
      return block_size_;
    }

    std::size_t NumColumns() const
    {
      return block_size_;
    }

    std::size_t FlatBlockSize() const
    {
      return number_of_non_zero_elements_per_block_;
    }

    /// @brief Set every matrix element to a given value
    /// @param val Value to set each element to
    void Fill(T val)
    {
      std::fill(data_.begin(), data_.end(), val);
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

    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix& matrix)
    {
      for (std::size_t i = 0; i < matrix.number_of_blocks_; ++i)
      {
        os << "Block " << i << std::endl;
        for (std::size_t j = 0; j < matrix.block_size_; ++j)
        {
          for (std::size_t k = 0; k < matrix.block_size_ - 1; ++k)
          {
            if (matrix.IsZero(j, k))
              os << "0,";
            else
              os << matrix[i][j][k] << ',';
          }
          if (matrix.IsZero(j, matrix.block_size_ - 1))
            os << "0" << std::endl;
          else
            os << matrix[i][j][matrix.block_size_ - 1] << std::endl;
        }
      }
      return os;
    }

    /// @brief Print the sparse matrix with row index, column index, and non-zero value; useful to test other linear algebra libraries
    /// @param os Output stream to print to, defaults to std::cout
    void PrintNonZeroElement(std::ostream& os) const
    {
      for (std::size_t i = 0; i < number_of_blocks_; ++i)
      {
        os << "Block " << i << std::endl;
        for (std::size_t j = 0; j < block_size_; ++j)
          for (std::size_t k = 0; k < block_size_; ++k)
            if (!this->IsZero(j, k)) os << j << ", " << k << ", " << (*this)[i][j][k] << std::endl;
          
      }
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

    SparseMatrixBuilder& SetNumberOfBlocks(std::size_t n)
    {
      number_of_blocks_ = n;
      return *this;
    }

    SparseMatrixBuilder& WithElement(std::size_t x, std::size_t y)
    {
      if (x >= block_size_ || y >= block_size_)
        throw std::system_error(make_error_code(MicmMatrixErrc::ElementOutOfRange));
      non_zero_elements_.insert(std::make_pair(x, y));
      return *this;
    }

    SparseMatrixBuilder& InitialValue(T inital_value)
    {
      initial_value_ = inital_value;
      return *this;
    }

    std::size_t NumberOfElements() const
    {
      return non_zero_elements_.size() * number_of_blocks_;
    }
  };

}  // namespace micm
