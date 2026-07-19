// Copyright (C) 2025-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/micm_exception.hpp>
#include <micm/util/types.hpp>
#include <micm/util/view_category.hpp>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iterator>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace micm
{

  /// @brief Defines the ordering of SparseMatrix object data in Compressed Sparse Column format
  ///
  /// Data is stored with blocks in the block diagonal matrix as the highest
  /// level structure, then by column, then by non-zero rows in each column.
  class SparseMatrixStandardOrderingCompressedSparseColumn
  {
    /// @brief Column ids of each non-zero element in a block
    std::vector<Index> column_ids_;
    /// @brief Start and end indices of each column in a block in column_ids_
    std::vector<Index> column_start_;
    /// @brief Indices along the diagonal of each block
    std::vector<Index> diagonal_ids_;

   protected:
    SparseMatrixStandardOrderingCompressedSparseColumn() = default;

    SparseMatrixStandardOrderingCompressedSparseColumn(
        Index number_of_blocks,
        Index block_size,
        const std::set<std::pair<Index, Index>>& non_zero_elements)
        : column_ids_(ColumnIdsVector(block_size, non_zero_elements)),
          column_start_(ColumnStartVector(block_size, non_zero_elements)),
          diagonal_ids_(DiagonalIndices(number_of_blocks, 0))
    {
    }

    SparseMatrixStandardOrderingCompressedSparseColumn& operator=(
        const std::tuple<Index, Index, std::set<std::pair<Index, Index>>>& number_size_elements)
    {
      column_ids_ = ColumnIdsVector(std::get<1>(number_size_elements), std::get<2>(number_size_elements));
      column_start_ = ColumnStartVector(std::get<1>(number_size_elements), std::get<2>(number_size_elements));
      diagonal_ids_ = DiagonalIndices(std::get<0>(number_size_elements), 0);
      return *this;
    }

    /// @brief Returns the size of the compressed data vector
    /// @param number_of_blocks Number of block sub-matrices in the overall matrix
    /// @return Size of the compressed data vector
    Index VectorSize(Index number_of_blocks) const
    {
      return number_of_blocks * column_ids_.size();
    }

    /// @brief Returns the index for a particular element in the compressed data vector
    /// @param number_of_blocks Total number of block sub-matrices in the overall matrix
    /// @param block Index of the block sub-matrix
    /// @param row Index of the row in the block
    /// @param column Index of the column in the block
    /// @return Index of the element in the compressed data vector
    Index VectorIndex(Index number_of_blocks, Index block, Index row, Index column) const
    {
      if (column >= column_start_.size() - 1 || row >= column_start_.size() - 1 || block >= number_of_blocks)
      {
        throw MicmException(MICM_ERROR_CATEGORY_MATRIX, MICM_MATRIX_ERROR_CODE_ELEMENT_OUT_OF_RANGE, "Element out of range");
      }
      auto begin = std::next(column_ids_.begin(), column_start_[column]);
      auto end = std::next(column_ids_.begin(), column_start_[column + 1]);
      auto elem = std::find(begin, end, row);
      if (elem == end)
      {
        throw MicmException(MICM_ERROR_CATEGORY_MATRIX, MICM_MATRIX_ERROR_CODE_ZERO_ELEMENT_ACCESS, "Zero element access");
      }
      return Index{ (elem - column_ids_.begin()) + block * column_ids_.size() };
    }

    /// @brief Adds a value to the diagonal of each block in the compressed data vector
    /// @param number_of_blocks Total number of block sub-matrices in the overall matrix
    /// @param data Compressed data vector
    /// @param value Value to add to the diagonal
    void AddToDiagonal(const Index number_of_blocks, auto& data, auto value) const
    {
      for (Index block_start = 0; block_start < number_of_blocks * column_ids_.size();
           block_start += column_ids_.size())
      {
        for (const auto& i : diagonal_ids_)
        {
          data[block_start + i] += value;
        }
      }
    }

    /// @brief Convert row and column indices to vector index within a block
    /// @param row The row index
    /// @param col The column index
    /// @return The index of the nth non-zero element within a block (0-based)
    Index VectorIndexFromRowColumn(Index row, Index col) const
    {
      if (col >= column_start_.size() - 1 || row >= column_start_.size() - 1)
      {
        throw MicmException(MICM_ERROR_CATEGORY_MATRIX, MICM_MATRIX_ERROR_CODE_ELEMENT_OUT_OF_RANGE, "Element out of range");
      }
      auto begin = std::next(column_ids_.begin(), column_start_[col]);
      auto end = std::next(column_ids_.begin(), column_start_[col + 1]);
      auto elem = std::find(begin, end, row);
      if (elem == end)
      {
        throw MicmException(MICM_ERROR_CATEGORY_MATRIX, MICM_MATRIX_ERROR_CODE_ZERO_ELEMENT_ACCESS, "Zero element access");
      }
      return std::distance(column_ids_.begin(), elem);
    }

    /// @brief Extract element position from VectorIndex(0, row, col) result
    /// @param vector_index_block_zero The result of VectorIndex(0, row, col)
    /// @return The element position (0 to number_of_non_zero_elements-1)
    Index ElementPositionFromVectorIndex(Index vector_index_block_zero) const
    {
      // For standard ordering: VectorIndex(0, row, col) = elem_position
      return vector_index_block_zero;
    }

    /// @brief Returns the indices along the diagonal of each block
    /// @param number_of_blocks Number of block sub-matrices in the overall matrix
    /// @param block_id Block index
    /// @return Vector of indices of non-zero diagonal elements
    std::vector<Index> DiagonalIndices(const Index number_of_blocks, const Index block_id) const
    {
      std::vector<Index> indices;
      indices.reserve(column_start_.size() - 1);
      for (Index i = 0; i < column_start_.size() - 1; ++i)
      {
        if (!IsZero(i, i))
        {
          indices.push_back(VectorIndex(number_of_blocks, block_id, i, i));
        }
      }
      return indices;
    }

   public:
    /// @brief A block-local temporary variable with its own storage
    /// For standard ordering: single value
    template<typename T>
    class BlockVariable
    {
     public:
      using category = BlockVariableTag;
      BlockVariable() = default;

      T& Get()
      {
        return storage_;
      }
      const T& Get() const
      {
        return storage_;
      }

     private:
      T storage_;
    };

    /// @brief ConstGroupView provides a const view of a single group of blocks for iteration
    /// For standard ordering: L=1, so each group contains 1 block
    template<typename SparseMatrixType>
    class ConstGroupView
    {
     private:
      const SparseMatrixType& matrix_;
      Index group_;

      /// @brief Get element from sparse matrix ConstBlockView
      template<SparseMatrixBlockView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg) const
      {
        auto* source_matrix = arg.GetMatrix();
        Index elem_position = arg.ElementPosition();
        Index num_non_zero = source_matrix->FlatBlockSize();
        // For standard ordering: block * num_non_zero + elem_position
        Index data_index = group_ * num_non_zero + elem_position;
        return source_matrix->AsVector()[data_index];
      }

      /// @brief Get element from Matrix or VectorMatrix ConstColumnView
      /// For standard ordering: compatible with standard Matrix or VectorMatrix with L=1
      template<DenseMatrixColumnView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg) const
      {
        auto* source_matrix = arg.GetMatrix();
        // Verify L=1 for VectorMatrix types
        if constexpr (requires { source_matrix->GroupVectorSize(); })
        {
          assert(
              source_matrix->GroupVectorSize() == 1 &&
              "Standard ordering sparse matrices require L=1 (use VectorMatrix<1>)");
        }
        return source_matrix->AsVector()[group_ * source_matrix->NumColumns() + arg.ColumnIndex()];
      }

      /// @brief Get element from BlockVariable
      template<BlockVariableView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg) const
      {
        return arg.Get();
      }

      /// @brief Get element from Vector-like
      template<VectorLike Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg) const
      {
        return arg[group_];
      }

     public:
      ConstGroupView(const SparseMatrixType& matrix, Index group)
          : matrix_(matrix),
            group_(group)
      {
      }

      auto GetConstBlockView(Index vector_index) const
      {
        return matrix_.GetConstBlockView(vector_index);
      }

      auto GetConstBlockView(Index row, Index col) const
      {
        return matrix_.GetConstBlockView(row, col);
      }

      auto GetBlockVariable() const
      {
        using T = typename SparseMatrixType::value_type;
        return BlockVariable<T>();
      }

      template<typename Func, typename... Args>
      void ForEachBlock(Func&& func, Args&&... args) const
      {
        // For standard ordering, only one block per group
        func(GetBlockElement(0, std::forward<Args>(args))...);
      }

      Index NumberOfBlocks() const
      {
        return matrix_.NumberOfBlocks();
      }
      Index NumRows() const
      {
        return matrix_.NumRows();
      }
      Index NumColumns() const
      {
        return matrix_.NumColumns();
      }
    };

    /// @brief GroupView provides a view of a single group of blocks for iteration
    /// For standard ordering: L=1, so each group contains 1 block
    template<typename SparseMatrixType>
    class GroupView
    {
     private:
      SparseMatrixType& matrix_;
      Index group_;

      /// @brief Get element from sparse matrix BlockView
      template<SparseMatrixBlockView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg)
      {
        auto* source_matrix = arg.GetMatrix();
        Index elem_position = arg.ElementPosition();
        Index num_non_zero = source_matrix->FlatBlockSize();
        // For standard ordering: block * num_non_zero + elem_position
        Index data_index = group_ * num_non_zero + elem_position;
        return source_matrix->AsVector()[data_index];
      }

      /// @brief Get element from Matrix or VectorMatrix ColumnView
      /// For standard ordering: compatible with standard Matrix or VectorMatrix with L=1
      template<DenseMatrixColumnView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg)
      {
        auto* source_matrix = arg.GetMatrix();
        // Verify L=1 for VectorMatrix types
        if constexpr (requires { source_matrix->GroupVectorSize(); })
        {
          assert(
              source_matrix->GroupVectorSize() == 1 &&
              "Standard ordering sparse matrices require L=1 (use VectorMatrix<1>)");
        }
        return source_matrix->AsVector()[group_ * source_matrix->NumColumns() + arg.ColumnIndex()];
      }

      /// @brief Get element from BlockVariable
      template<BlockVariableView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg)
      {
        return arg.Get();
      }

      /// @brief Get element from Vector-like
      template<VectorLike Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg)
      {
        return arg[group_];
      }

     public:
      GroupView(SparseMatrixType& matrix, Index group)
          : matrix_(matrix),
            group_(group)
      {
      }

      auto GetConstBlockView(Index vector_index) const
      {
        return matrix_.GetConstBlockView(vector_index);
      }

      auto GetConstBlockView(Index row, Index col) const
      {
        return matrix_.GetConstBlockView(row, col);
      }

      auto GetBlockView(Index vector_index)
      {
        return matrix_.GetBlockView(vector_index);
      }

      auto GetBlockView(Index row, Index col)
      {
        return matrix_.GetBlockView(row, col);
      }

      auto GetBlockVariable()
      {
        using T = typename SparseMatrixType::value_type;
        return BlockVariable<T>();
      }

      template<typename Func, typename... Args>
      void ForEachBlock(Func&& func, Args&&... args)
      {
        // For standard ordering, only one block per group
        func(GetBlockElement(0, std::forward<Args>(args))...);
      }

      Index NumberOfBlocks() const
      {
        return matrix_.NumberOfBlocks();
      }
      Index NumRows() const
      {
        return matrix_.NumRows();
      }
      Index NumColumns() const
      {
        return matrix_.NumColumns();
      }
    };

    /// @brief Returns the number of blocks included in each group of blocks
    /// @return Number of blocks in each group (1 for standard ordering)
    static constexpr Index GroupVectorSize()
    {
      return 1;
    }

    /// @brief Returns the size of each group of blocks in the compressed data vector
    /// @return Size of each group of blocks
    Index GroupSize(Index number_of_non_zero_elements) const
    {
      return number_of_non_zero_elements;
    }

    /// @brief Returns the total number of groups of blocks in the compressed data
    /// @param number_of_blocks Total number of block sub-matrices in the overall matrix
    /// @return Number of groups of blocks (equal to number_of_blocks for standard ordering)
    Index NumberOfGroups(Index number_of_blocks) const
    {
      return number_of_blocks;
    }

   private:
    /// @brief Returns the column ids of each non-zero element in a block
    /// @param block_size Number of rows or columns for each block
    /// @param non_zero_elements Set of non-zero elements in the matrix
    static std::vector<Index> ColumnIdsVector(
        const Index block_size,
        const std::set<std::pair<Index, Index>>& non_zero_elements)
    {
      std::vector<Index> ids;
      ids.reserve(non_zero_elements.size());
      std::set<std::pair<Index, Index>> column_ordered_pairs;
      for (const auto& elem : non_zero_elements)
      {
        column_ordered_pairs.insert(std::make_pair(elem.second, elem.first));
      }
      std::transform(
          column_ordered_pairs.begin(),
          column_ordered_pairs.end(),
          std::back_inserter(ids),
          [](const std::pair<Index, Index>& elem) { return elem.second; });
      return ids;
    }

    /// @brief Returns the start and end indices of each column in a block in column_ids_
    /// @param block_size Number of rows or columns for each block
    /// @param non_zero_elements Set of non-zero elements in the matrix
    static std::vector<Index> ColumnStartVector(
        const Index block_size,
        const std::set<std::pair<Index, Index>>& non_zero_elements)
    {
      std::vector<Index> starts(block_size + 1, 0);
      Index total_elem = 0;
      Index curr_col = 0;
      std::set<std::pair<Index, Index>> column_ordered_pairs;
      for (const auto& elem : non_zero_elements)
      {
        column_ordered_pairs.insert(std::make_pair(elem.second, elem.first));
      }
      for (const auto& elem : column_ordered_pairs)
      {
        while (curr_col < elem.first)
        {
          starts[(curr_col++) + 1] = total_elem;
        }
        ++total_elem;
      }
      // Fill all remaining entries from curr_col + 1 to block_size
      for (Index i = curr_col + 1; i <= block_size; ++i)
      {
        starts[i] = total_elem;
      }
      return starts;
    }

   public:
    /// @brief Returns whether a particular element is always zero
    /// @param row Row index
    /// @param column Column index
    /// @return true if the element is always zero, false otherwise
    bool IsZero(Index row, Index column) const
    {
      if (column >= column_start_.size() - 1 || row >= column_start_.size() - 1)
      {
        throw MicmException(MICM_ERROR_CATEGORY_MATRIX, MICM_MATRIX_ERROR_CODE_ELEMENT_OUT_OF_RANGE, "Element out of range");
      }
      auto begin = std::next(column_ids_.begin(), column_start_[column]);
      auto end = std::next(column_ids_.begin(), column_start_[column + 1]);
      auto elem = std::find(begin, end, row);
      return (elem == end);
    }
  };

}  // namespace micm