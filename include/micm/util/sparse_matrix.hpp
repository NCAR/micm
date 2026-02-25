// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/matrix_error.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/view_category.hpp>

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
      using category = SparseMatrixBlockViewTag;
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
      using category = SparseMatrixBlockViewTag;
      std::size_t RowIndex() const { return row_index_; }
      std::size_t ColumnIndex() const { return column_index_; }
      SparseMatrix* GetMatrix() { return matrix_; }
    };

    /// @brief Alias for the ordering policy's BlockVariable type
    using BlockVariable = typename OrderingPolicy::template BlockVariable<T>;

    /// @brief Alias for the ordering policy's ConstGroupView type
    using ConstGroupView = typename OrderingPolicy::template ConstGroupView<SparseMatrix>;

    /// @brief Alias for the ordering policy's GroupView type
    using GroupView = typename OrderingPolicy::template GroupView<SparseMatrix>;

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

    /// @brief Create a function that can be applied to sparse matrices and vectors
    /// 
    /// Creates a reusable callable that validates dimensions and applies a user function
    /// across block groups. The function iterates over groups of L blocks at a time, where L
    /// is determined by the OrderingPolicy::GroupVectorSize(). Supports mixing sparse matrices,
    /// dense matrices, and vector-like types.
    /// 
    /// @tparam Func The lambda/function type
    /// @tparam Args The matrix and vector types (can mix SparseMatrix, VectorMatrix, Matrix, and vectors)
    /// @param func The function to wrap - receives GroupView objects for matrices and forwarded vectors
    /// @param args The matrices and vectors to validate and capture dimensions from
    /// @return A callable that validates dimensions and applies the function
    /// 
    /// @note Validation occurs in two phases:
    ///   1. At function creation: Validates matrix dimensions, vector sizes, and ordering compatibility
    ///   2. At invocation: Re-validates dimensions in case matrices/vectors were resized
    /// 
    /// @note Column/Block view creation happens inside user lambda and is validated
    ///       at invocation time, not at function creation time. Ensure all view indices
    ///       are within matrix bounds to avoid runtime errors.
    /// 
    /// @throws std::system_error if matrices have incompatible orderings (different L values),
    ///         mismatched block counts, or vectors have wrong sizes
    template<typename Func, typename... Args>
    static auto Function(Func&& func, Args&... args)
    {
      // Validate that all matrices have compatible ordering (same L value)
      // Get L from this sparse matrix's ordering policy
      constexpr std::size_t expected_L = OrderingPolicy::GroupVectorSize();
      
      
      // Check each argument: matrices must have compatible L, vectors are skipped
      std::size_t index = 0;
      ([&](auto& arg) {
        using ArgType = std::remove_cvref_t<decltype(arg)>;
        
        // Only check L for matrix types (not vectors)
        if constexpr (!VectorLike<ArgType>)
        {
          // Get the L value for this matrix using the type trait
          constexpr std::size_t matrix_L = GroupVectorSize_v<std::decay_t<decltype(arg)>>;
          
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
        }
        ++index;
      }(args), ...);
      
      // Validate that all matrices have compatible dimensions and vectors have matching sizes
      // For sparse matrices: use NumberOfBlocks()
      // For dense matrices: use NumRows() (blocks correspond to rows)
      // For vectors: use size() (should match block count)
      std::size_t num_blocks = 0;
      std::array<std::size_t, sizeof...(Args)> block_sizes_or_vec_sizes{};
      std::array<bool, sizeof...(Args)> is_sparse{};
      std::array<bool, sizeof...(Args)> is_vector{};
      bool found_first = false;
      index = 0;
      
      ([&](auto& arg) {
        using ArgType = std::remove_cvref_t<decltype(arg)>;
        
        // Check if this is a vector-like type
        constexpr bool is_vector_type = VectorLike<ArgType>;
        is_vector[index] = is_vector_type;
        
        if constexpr (is_vector_type)
        {
          // This is a vector
          is_sparse[index] = false;
          if (!found_first)
          {
            num_blocks = arg.size();
            found_first = true;
          }
          else if (arg.size() != num_blocks)
          {
            throw std::system_error(
                make_error_code(MicmMatrixErrc::InvalidVector),
                "Vector size " + std::to_string(arg.size()) + 
                    " does not match expected block count " + std::to_string(num_blocks));
          }
          block_sizes_or_vec_sizes[index] = arg.size();
        }
        else
        {
          // This is a matrix
          constexpr bool is_sparse_matrix = SparseMatrixConcept<ArgType>;
          is_sparse[index] = is_sparse_matrix;
          
          if (!found_first)
          {
            if constexpr (is_sparse_matrix)
            {
              num_blocks = arg.NumberOfBlocks();
            }
            else
            {
              // For dense matrices, rows correspond to blocks
              num_blocks = arg.NumRows();
            }
            found_first = true;
          }
          else
          {
            std::size_t arg_blocks;
            if constexpr (is_sparse_matrix)
            {
              arg_blocks = arg.NumberOfBlocks();
            }
            else
            {
              // For dense matrices, rows correspond to blocks
              arg_blocks = arg.NumRows();
            }
            
            if (arg_blocks != num_blocks)
            {
              throw std::system_error(
                  make_error_code(MicmMatrixErrc::InvalidVector),
                  "All matrices must have the same number of blocks/rows. Expected " + std::to_string(num_blocks) +
                      " but got " + std::to_string(arg_blocks));
            }
          }
          
          // Store block size (for sparse) or number of columns (for dense)
          block_sizes_or_vec_sizes[index] = arg.NumRows();  // For sparse: block size; for dense: also NumRows
        }
        ++index;
      }(args), ...);

      // Return a callable that validates dimensions on invocation and applies the function
      return [func = std::forward<Func>(func), num_blocks, block_sizes_or_vec_sizes, is_sparse, is_vector](Args&... invoked_args) {
        std::size_t idx = 0;
        ([&](auto& arg) {
          using ArgType = std::remove_cvref_t<decltype(arg)>;
          
          if constexpr (VectorLike<ArgType>)
          {
            // Validate vector size
            if (arg.size() != block_sizes_or_vec_sizes[idx])
            {
              throw std::system_error(
                  make_error_code(MicmMatrixErrc::InvalidVector),
                  "Vector dimensions do not match. Expected " + std::to_string(block_sizes_or_vec_sizes[idx]) + 
                      " elements but got " + std::to_string(arg.size()));
            }
          }
          else
          {
            // Validate matrix dimensions
            std::size_t arg_blocks;
            constexpr bool is_sparse_matrix = SparseMatrixConcept<ArgType>;
            
            if constexpr (is_sparse_matrix)
            {
              arg_blocks = arg.NumberOfBlocks();
            }
            else
            {
              // For dense matrices, rows correspond to blocks
              arg_blocks = arg.NumRows();
            }
            
            if (arg_blocks != num_blocks)
            {
              throw std::system_error(
                  make_error_code(MicmMatrixErrc::InvalidVector),
                  "Matrix dimensions do not match. Expected " + std::to_string(num_blocks) + 
                      " blocks/rows but got " + std::to_string(arg_blocks));
            }
            if (arg.NumRows() != block_sizes_or_vec_sizes[idx])
            {
              throw std::system_error(
                  make_error_code(MicmMatrixErrc::InvalidVector),
                  "Matrix block size/rows does not match. Expected " + std::to_string(block_sizes_or_vec_sizes[idx]) + 
                      " but got " + std::to_string(arg.NumRows()));
            }
          }
          ++idx;
        }(invoked_args), ...);
        
        // Get the group vector size from the OrderingPolicy (compile-time constant)
        // For standard ordering: L = 1
        // For vector ordering: L > 1
        constexpr std::size_t L = OrderingPolicy::GroupVectorSize();
        
        // Iterate over groups, processing L blocks at a time
        std::size_t num_groups = (num_blocks + L - 1) / L;  // Ceiling division
        for (std::size_t group = 0; group < num_groups; ++group)
        {
          // For matrices: use ConstGroupView if const, otherwise GroupView
          // For vectors: forward them directly
          func([&]() -> decltype(auto) {
            using ArgType = std::remove_reference_t<decltype(invoked_args)>;
            if constexpr (VectorLike<std::remove_cvref_t<ArgType>>)
            {
              // Vector: just forward it
              return std::forward<decltype(invoked_args)>(invoked_args);
            }
            else
            {
              // Matrix: create appropriate GroupView
              if constexpr (std::is_const_v<ArgType>)
              {
                return typename std::decay_t<Args>::ConstGroupView(invoked_args, group);
              }
              else
              {
                return typename std::decay_t<Args>::GroupView(invoked_args, group);
              }
            }
          }()...);
        }
      };
    }

   private:
    /// @brief Get an element reference for a block (BlockView)
    template<SparseMatrixBlockView Arg>
    [[gnu::always_inline]]
    inline decltype(auto) GetBlockElement(std::size_t block, Arg&& arg)
    {
      auto* source_matrix = arg.GetMatrix();
      return source_matrix->data_[source_matrix->VectorIndex(block, arg.RowIndex(), arg.ColumnIndex())];
    }

    /// @brief Get an element reference for a block (BlockVariable)
    template<BlockVariableView Arg>
    [[gnu::always_inline]]
    inline decltype(auto) GetBlockElement(std::size_t block, Arg&& arg)
    {
      return arg.Get();
    }

    /// @brief Get an element reference for a block (Vector-like)
    template<VectorLike Arg>
    [[gnu::always_inline]]
    inline decltype(auto) GetBlockElement(std::size_t block, Arg&& arg)
    {
      return arg[block];
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
