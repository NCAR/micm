// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/matrix_error.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
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

  /// @brief Type trait to extract GroupVectorSize (L) from matrix types at compile-time
  /// Default: L=1 for types without GroupVectorSize
  template<typename T>
  struct GroupVectorSize : std::integral_constant<std::size_t, 1> {};

  /// Specialization for types with static GroupVectorSize method
  template<typename T>
    requires requires { T::GroupVectorSize(); }
  struct GroupVectorSize<T> : std::integral_constant<std::size_t, T::GroupVectorSize()> {};

  /// Helper variable template
  template<typename T>
  inline constexpr std::size_t GroupVectorSize_v = GroupVectorSize<T>::value;

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

    /// @brief A lightweight descriptor for a const block element in a sparse matrix
    class ConstBlockView
    {
      friend class SparseMatrix;
      const SparseMatrix* matrix_;
      std::size_t row_index_;
      std::size_t column_index_;
      
      explicit ConstBlockView(const SparseMatrix* matrix, std::size_t row_index, std::size_t column_index)
          : matrix_(matrix),
            row_index_(row_index),
            column_index_(column_index)
      {
      }

     public:
      std::size_t RowIndex() const { return row_index_; }
      std::size_t ColumnIndex() const { return column_index_; }
      const SparseMatrix* GetMatrix() const { return matrix_; }
    };

    /// @brief A lightweight descriptor for a mutable block element in a sparse matrix
    class BlockView
    {
      friend class SparseMatrix;
      SparseMatrix* matrix_;
      std::size_t row_index_;
      std::size_t column_index_;
      
      explicit BlockView(SparseMatrix* matrix, std::size_t row_index, std::size_t column_index)
          : matrix_(matrix),
            row_index_(row_index),
            column_index_(column_index)
      {
      }

     public:
      std::size_t RowIndex() const { return row_index_; }
      std::size_t ColumnIndex() const { return column_index_; }
      SparseMatrix* GetMatrix() { return matrix_; }
    };

    /// @brief A block-local temporary variable with its own storage
    /// For standard ordering: single value
    /// For vector ordering: array of L values
    class BlockVariable
    {
      friend class SparseMatrix;
      friend class GroupView;
      
      // Use conditional storage based on ordering policy
      static constexpr std::size_t L = OrderingPolicy::GroupVectorSize();
      typename std::conditional<(L > 1), std::array<T, L>, T>::type storage_;
      
     public:
      BlockVariable() = default;
      
      auto& Get() { return storage_; }
      const auto& Get() const { return storage_; }
    };

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
    template<typename> friend class Matrix;

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

    /// @brief Constructs a SparseMatrix from a given builder and optional indexing mode.
    ///        Initializes the SparseMatrix using the provided SparseMatrixBuilder, which defines
    ///        the matrix structure, block size, and non-zero elements. Optionally, the constructor
    ///        can be used in "indexing only" mode, where the data storage is not allocated.
    /// @tparam T The type of the matrix elements.
    /// @tparam OrderingPolicy The policy class that defines the ordering and storage of elements.
    /// @param builder The builder object containing matrix configuration and initial values.
    /// @param indexing_only If true, only indexing structures are initialized and data storage is omitted.
    SparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder, bool indexing_only = false)
        : OrderingPolicy(builder.number_of_blocks_, builder.block_size_, builder.non_zero_elements_),
          number_of_blocks_(builder.number_of_blocks_),
          block_size_(builder.block_size_),
          number_of_non_zero_elements_per_block_(builder.non_zero_elements_.size()),
          data_((indexing_only ? 0 : OrderingPolicy::VectorSize(number_of_blocks_)), builder.initial_value_)
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

    /// @brief Print the sparse matrix with row index, column index, and non-zero value; useful to test other linear algebra
    /// libraries
    /// @param os Output stream to print to, defaults to std::cout
    void PrintNonZeroElements(std::ostream& os) const
    {
      for (std::size_t i = 0; i < number_of_blocks_; ++i)
      {
        os << "Block " << i << std::endl;
        for (std::size_t j = 0; j < block_size_; ++j)
          for (std::size_t k = 0; k < block_size_; ++k)
            if (!this->IsZero(j, k))
              os << j << ", " << k << ", " << (*this)[i][j][k] << std::endl;
      }
    }

    /// @brief Create a const block view for accessing a block element
    /// @param row The row index of the block element
    /// @param col The column index of the block element
    /// @return A ConstBlockView descriptor
    ConstBlockView GetConstBlockView(std::size_t row, std::size_t col) const
    {
      if (row >= block_size_ || col >= block_size_)
      {
        throw std::system_error(
            make_error_code(MicmMatrixErrc::ElementOutOfRange),
            "Block element (" + std::to_string(row) + "," + std::to_string(col) + 
            ") out of range for matrix with block size " + std::to_string(block_size_));
      }
      if (this->IsZero(row, col))
      {
        throw std::system_error(
            make_error_code(MicmMatrixErrc::ZeroElementAccess),
            "Cannot create view for zero block element (" + std::to_string(row) + "," + 
            std::to_string(col) + ")");
      }
      return ConstBlockView(this, row, col);
    }

    /// @brief Create a mutable block view for accessing a block element
    /// @param row The row index of the block element
    /// @param col The column index of the block element
    /// @return A BlockView descriptor
    BlockView GetBlockView(std::size_t row, std::size_t col)
    {
      if (row >= block_size_ || col >= block_size_)
      {
        throw std::system_error(
            make_error_code(MicmMatrixErrc::ElementOutOfRange),
            "Block element (" + std::to_string(row) + "," + std::to_string(col) + 
            ") out of range for matrix with block size " + std::to_string(block_size_));
      }
      if (this->IsZero(row, col))
      {
        throw std::system_error(
            make_error_code(MicmMatrixErrc::ZeroElementAccess),
            "Cannot create view for zero block element (" + std::to_string(row) + "," + 
            std::to_string(col) + ")");
      }
      return BlockView(this, row, col);
    }

    /// @brief Get a block variable with persistent storage for temporary values
    /// @return A BlockVariable with stack-allocated storage
    BlockVariable GetBlockVariable()
    {
      return BlockVariable();
    }

    /// @brief Apply a function to each block of the matrix
    /// @tparam Func The lambda/function type
    /// @tparam Args The types of the block view arguments
    /// @param func The function to apply to each block
    /// @param args Block views or block variables
    template<typename Func, typename... Args>
    void ForEachBlock(Func&& func, Args&&... args)
    {
      for (std::size_t block = 0; block < number_of_blocks_; ++block)
      {
        func(GetBlockElement(block, args)...);
      }
    }

    /// @brief ConstGroupView provides a const view of a single group of blocks for iteration
    class ConstGroupView
    {
     private:
      const SparseMatrix& matrix_;
      std::size_t group_;
      std::size_t num_blocks_in_group_;  // May be < L for the last group in vector ordering

      /// @brief Get a const element reference for a specific block in this group
      template<typename Arg>
      [[gnu::always_inline]]
      inline decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg) const
      {
        // Calculate the actual block index from group and block_in_group
        std::size_t block = group_ * OrderingPolicy::GroupVectorSize() + block_in_group;
        
        // Check if Arg has RowIndex() method (ConstBlockView from sparse matrix)
        if constexpr (requires { arg.RowIndex(); arg.ColumnIndex(); })
        {
          // It's a ConstBlockView type from a sparse matrix
          auto* source_matrix = arg.GetMatrix();
          return source_matrix->data_[source_matrix->VectorIndex(block, arg.RowIndex(), arg.ColumnIndex())];
        }
        // Check if Arg has ColumnIndex() method but not RowIndex() (ConstColumnView from dense matrix)
        else if constexpr (requires { arg.ColumnIndex(); } && !requires { arg.RowIndex(); })
        {
          // It's a ConstColumnView type from a dense matrix
          auto* source_matrix = arg.GetMatrix();
          
          // Check if this is a VectorMatrix (has GroupVectorSize method)
          if constexpr (requires { source_matrix->GroupVectorSize(); })
          {
            // VectorMatrix layout: data_[(group * y_dim + column) * L + row_in_group]
            std::size_t L = source_matrix->GroupVectorSize();
            std::size_t row = block;
            std::size_t row_group = row / L;
            std::size_t row_in_group = row % L;
            return source_matrix->data_[(row_group * source_matrix->NumColumns() + arg.ColumnIndex()) * L + row_in_group];
          }
          else
          {
            // Standard Matrix layout: data_[row * num_cols + col]
            return source_matrix->data_[block * source_matrix->NumColumns() + arg.ColumnIndex()];
          }
        }
        else if constexpr (requires { arg.Get(); })
        {
          // It's a BlockVariable, access the array element for vector ordering
          if constexpr (OrderingPolicy::GroupVectorSize() > 1)
          {
            // Vector ordering: BlockVariable has array storage
            return arg.Get()[block_in_group];
          }
          else
          {
            // Standard ordering: BlockVariable has single value storage
            return arg.Get();
          }
        }
        else
        {
          // Unknown type, just return it
          return arg;
        }
      }

     public:
      ConstGroupView(const SparseMatrix& matrix, std::size_t group)
          : matrix_(matrix), group_(group)
      {
        // Calculate how many blocks are in this group (compile-time L, runtime calculation)
        constexpr std::size_t L = OrderingPolicy::GroupVectorSize();
        std::size_t total_blocks = matrix_.NumberOfBlocks();
        std::size_t start_block = group * L;
        num_blocks_in_group_ = std::min(L, total_blocks - start_block);
      }

      auto GetConstBlockView(std::size_t row, std::size_t col) const
      {
        return matrix_.GetConstBlockView(row, col);
      }

      BlockVariable GetBlockVariable() const
      {
        return BlockVariable();
      }

      template<typename Func, typename... Args>
      void ForEachBlock(Func&& func, Args&&... args) const
      {
        // Tight loop over blocks in this group for vectorization
        for (std::size_t block_in_group = 0; block_in_group < num_blocks_in_group_; ++block_in_group)
        {
          func(GetBlockElement(block_in_group, std::forward<Args>(args))...);
        }
      }

      std::size_t NumberOfBlocks() const { return matrix_.NumberOfBlocks(); }
      std::size_t NumRows() const { return matrix_.NumRows(); }
      std::size_t NumColumns() const { return matrix_.NumColumns(); }
    };

    /// @brief GroupView provides a view of a single group of blocks for iteration
    class GroupView
    {
     private:
      SparseMatrix& matrix_;
      std::size_t group_;
      std::size_t num_blocks_in_group_;  // May be < L for the last group in vector ordering

      /// @brief Get an element reference for a specific block in this group
      template<typename Arg>
      [[gnu::always_inline]]
      inline decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg)
      {
        // Calculate the actual block index from group and block_in_group
        std::size_t block = group_ * OrderingPolicy::GroupVectorSize() + block_in_group;
        
        // Check if Arg has RowIndex() method (BlockView from sparse matrix)
        if constexpr (requires { arg.RowIndex(); arg.ColumnIndex(); })
        {
          // It's a BlockView type from a sparse matrix
          auto* source_matrix = arg.GetMatrix();
          return source_matrix->data_[source_matrix->VectorIndex(block, arg.RowIndex(), arg.ColumnIndex())];
        }
        // Check if Arg has ColumnIndex() method but not RowIndex() (ColumnView from dense matrix)
        else if constexpr (requires { arg.ColumnIndex(); } && !requires { arg.RowIndex(); })
        {
          // It's a ColumnView type from a dense matrix
          auto* source_matrix = arg.GetMatrix();
          
          // Check if this is a VectorMatrix (has GroupVectorSize method)
          if constexpr (requires { source_matrix->GroupVectorSize(); })
          {
            // VectorMatrix layout: data_[(group * y_dim + column) * L + row_in_group]
            std::size_t L = source_matrix->GroupVectorSize();
            std::size_t row = block;
            std::size_t row_group = row / L;
            std::size_t row_in_group = row % L;
            return source_matrix->data_[(row_group * source_matrix->NumColumns() + arg.ColumnIndex()) * L + row_in_group];
          }
          else
          {
            // Standard Matrix layout: data_[row * num_cols + col]
            return source_matrix->data_[block * source_matrix->NumColumns() + arg.ColumnIndex()];
          }
        }
        else if constexpr (requires { arg.Get(); })
        {
          // It's a BlockVariable, access the array element for vector ordering
          if constexpr (OrderingPolicy::GroupVectorSize() > 1)
          {
            // Vector ordering: BlockVariable has array storage
            return arg.Get()[block_in_group];
          }
          else
          {
            // Standard ordering: BlockVariable has single value storage
            return arg.Get();
          }
        }
        else
        {
          // Unknown type, just return it
          return arg;
        }
      }

     public:
      GroupView(SparseMatrix& matrix, std::size_t group)
          : matrix_(matrix), group_(group)
      {
        // Calculate how many blocks are in this group (compile-time L, runtime calculation)
        constexpr std::size_t L = OrderingPolicy::GroupVectorSize();
        std::size_t total_blocks = matrix_.NumberOfBlocks();
        std::size_t start_block = group * L;
        num_blocks_in_group_ = std::min(L, total_blocks - start_block);
      }

      auto GetConstBlockView(std::size_t row, std::size_t col) const
      {
        return matrix_.GetConstBlockView(row, col);
      }

      auto GetBlockView(std::size_t row, std::size_t col)
      {
        return matrix_.GetBlockView(row, col);
      }

      BlockVariable GetBlockVariable()
      {
        return BlockVariable();
      }

      template<typename Func, typename... Args>
      void ForEachBlock(Func&& func, Args&&... args)
      {
        // Tight loop over blocks in this group for vectorization
        for (std::size_t block_in_group = 0; block_in_group < num_blocks_in_group_; ++block_in_group)
        {
          func(GetBlockElement(block_in_group, std::forward<Args>(args))...);
        }
      }

      std::size_t NumberOfBlocks() const { return matrix_.NumberOfBlocks(); }
      std::size_t NumRows() const { return matrix_.NumRows(); }
      std::size_t NumColumns() const { return matrix_.NumColumns(); }
    };

    /// @brief Create a function that can be applied to sparse matrices
    /// 
    /// Creates a reusable callable that validates matrix dimensions and applies a user function
    /// across block groups. The function iterates over groups of L blocks at a time, where L
    /// is determined by the OrderingPolicy::GroupVectorSize().
    /// 
    /// @tparam Func The lambda/function type
    /// @tparam Matrices The matrix types (can mix SparseMatrix, VectorMatrix, and Matrix)
    /// @param func The function to wrap - receives GroupView objects for each matrix
    /// @param matrices The matrices to validate and capture dimensions from
    /// @return A callable that validates dimensions and applies the function
    /// 
    /// @note Validation occurs in two phases:
    ///   1. At function creation: Validates matrix dimensions and ordering compatibility
    ///   2. At invocation: Re-validates dimensions in case matrices were resized
    /// 
    /// @note Column/Block view creation happens inside user lambda and is validated
    ///       at invocation time, not at function creation time. Ensure all view indices
    ///       are within matrix bounds to avoid runtime errors.
    /// 
    /// @throws std::system_error if matrices have incompatible orderings (different L values)
    ///         or mismatched block counts
    template<typename Func, typename... Matrices>
    static auto Function(Func&& func, Matrices&... matrices)
    {
      // Validate that all matrices have compatible ordering (same L value)
      // Get L from this sparse matrix's ordering policy
      constexpr std::size_t expected_L = OrderingPolicy::GroupVectorSize();
      
      // Check each matrix has compatible L
      std::size_t index = 0;
      ([&](auto& matrix) {
        // Get the L value for this matrix using the type trait
        constexpr std::size_t matrix_L = GroupVectorSize_v<std::decay_t<decltype(matrix)>>;
        
        if (matrix_L != expected_L)
        {
          throw std::system_error(
              make_error_code(MicmMatrixErrc::InvalidVector),
              "Incompatible matrix orderings: Matrix " + std::to_string(index) + 
              " has GroupVectorSize=" + std::to_string(matrix_L) +
              " but expected " + std::to_string(expected_L) +
              ". Cannot mix standard-ordered (L=1) and vector-ordered (L>1) matrices, " +
              "or vector-ordered matrices with different L values.");
        }
        ++index;
      }(matrices), ...);
      
      // Validate that all matrices have compatible dimensions
      // For sparse matrices: use NumberOfBlocks()
      // For dense matrices: use NumRows() (blocks correspond to rows)
      std::size_t num_blocks = 0;
      std::array<std::size_t, sizeof...(Matrices)> block_sizes{};
      std::array<bool, sizeof...(Matrices)> is_sparse{};
      index = 0;
      
      ([&](auto& matrix) {
        // Check if this matrix is sparse (has NumberOfBlocks method)
        constexpr bool has_number_of_blocks = requires { matrix.NumberOfBlocks(); };
        is_sparse[index] = has_number_of_blocks;
        
        if (index == 0)
        {
          if constexpr (has_number_of_blocks)
          {
            num_blocks = matrix.NumberOfBlocks();
          }
          else
          {
            // For dense matrices, rows correspond to blocks
            num_blocks = matrix.NumRows();
          }
        }
        else
        {
          std::size_t matrix_blocks;
          if constexpr (has_number_of_blocks)
          {
            matrix_blocks = matrix.NumberOfBlocks();
          }
          else
          {
            // For dense matrices, rows correspond to blocks
            matrix_blocks = matrix.NumRows();
          }
          
          if (matrix_blocks != num_blocks)
          {
            throw std::system_error(
                make_error_code(MicmMatrixErrc::InvalidVector),
                "All matrices must have the same number of blocks/rows. Expected " + std::to_string(num_blocks) +
                    " but got " + std::to_string(matrix_blocks));
          }
        }
        
        // Store block size (for sparse) or number of columns (for dense)
        block_sizes[index] = matrix.NumRows();  // For sparse: block size; for dense: also NumRows
        ++index;
      }(matrices), ...);

      // Return a callable that validates dimensions on invocation and applies the function
      return [func = std::forward<Func>(func), num_blocks, block_sizes, is_sparse](Matrices&... invoked_matrices) {
        std::size_t idx = 0;
        ([&](auto& matrix) {
          std::size_t matrix_blocks;
          constexpr bool has_number_of_blocks = requires { matrix.NumberOfBlocks(); };
          
          if constexpr (has_number_of_blocks)
          {
            matrix_blocks = matrix.NumberOfBlocks();
          }
          else
          {
            // For dense matrices, rows correspond to blocks
            matrix_blocks = matrix.NumRows();
          }
          
          if (matrix_blocks != num_blocks)
          {
            throw std::system_error(
                make_error_code(MicmMatrixErrc::InvalidVector),
                "Matrix dimensions do not match. Expected " + std::to_string(num_blocks) + 
                    " blocks/rows but got " + std::to_string(matrix_blocks));
          }
          if (matrix.NumRows() != block_sizes[idx])
          {
            throw std::system_error(
                make_error_code(MicmMatrixErrc::InvalidVector),
                "Matrix block size/rows does not match. Expected " + std::to_string(block_sizes[idx]) + 
                    " but got " + std::to_string(matrix.NumRows()));
          }
          ++idx;
        }(invoked_matrices), ...);
        
        // Get the group vector size from the OrderingPolicy (compile-time constant)
        // For standard ordering: L = 1
        // For vector ordering: L > 1
        constexpr std::size_t L = OrderingPolicy::GroupVectorSize();
        
        // Iterate over groups, processing L blocks at a time
        std::size_t num_groups = (num_blocks + L - 1) / L;  // Ceiling division
        for (std::size_t group = 0; group < num_groups; ++group)
        {
          // Use ConstGroupView if matrix is const, otherwise use GroupView
          func([&]() {
            using MatrixType = std::remove_reference_t<decltype(invoked_matrices)>;
            if constexpr (std::is_const_v<MatrixType>)
            {
              return typename std::decay_t<Matrices>::ConstGroupView(invoked_matrices, group);
            }
            else
            {
              return typename std::decay_t<Matrices>::GroupView(invoked_matrices, group);
            }
          }()...);
        }
      };
    }

   private:
    /// @brief Get an element reference for a block, handling BlockViews and BlockVariables
    template<typename Arg>
    [[gnu::always_inline]]
    inline decltype(auto) GetBlockElement(std::size_t block, Arg&& arg)
    {
      // Check if Arg has GetMatrix() method (BlockView from potentially different matrix)
      if constexpr (requires { arg.GetMatrix(); })
      {
        // It's a BlockView type, access the source matrix's data
        auto* source_matrix = arg.GetMatrix();
        return source_matrix->data_[source_matrix->VectorIndex(block, arg.RowIndex(), arg.ColumnIndex())];
      }
      else if constexpr (requires { arg.Get(); })
      {
        // It's a BlockVariable, return reference to the single storage value
        return arg.Get();
      }
      else
      {
        // Unknown type, just return it
        return arg;
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
