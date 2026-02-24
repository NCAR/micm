// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/matrix_error.hpp>

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace micm
{

  /// @brief Defines the ordering of SparseMatrix object data in Compressed Sparse Row format
  ///
  /// Data is stored with blocks in the block diagonal matrix as the highest
  /// level structure, then by row, then by non-zero columns in each row.
  class SparseMatrixStandardOrderingCompressedSparseRow
  {
    /// @brief Row ids of each non-zero element in a block
    std::vector<std::size_t> row_ids_;
    /// @brief Start and end indices of each row in a block in row_ids_
    std::vector<std::size_t> row_start_;
    /// @brief Indices along the diagonal of each block
    std::vector<std::size_t> diagonal_ids_;

   protected:
    SparseMatrixStandardOrderingCompressedSparseRow() = default;

    SparseMatrixStandardOrderingCompressedSparseRow(
        std::size_t number_of_blocks,
        std::size_t block_size,
        std::set<std::pair<std::size_t, std::size_t>> non_zero_elements)
        : row_ids_(RowIdsVector(block_size, non_zero_elements)),
          row_start_(RowStartVector(block_size, non_zero_elements)),
          diagonal_ids_(DiagonalIndices(number_of_blocks, 0))
    {
    }

    SparseMatrixStandardOrderingCompressedSparseRow& operator=(
        const std::tuple<std::size_t, std::size_t, std::set<std::pair<std::size_t, std::size_t>>> number_size_elements)
    {
      row_ids_ = RowIdsVector(std::get<1>(number_size_elements), std::get<2>(number_size_elements));
      row_start_ = RowStartVector(std::get<1>(number_size_elements), std::get<2>(number_size_elements));
      diagonal_ids_ = DiagonalIndices(std::get<0>(number_size_elements), 0);
      return *this;
    }

    /// @brief Returns the size of the compressed data vector
    /// @param number_of_blocks Number of block sub-matrices in the overall matrix
    /// @return Size of the compressed data vector
    std::size_t VectorSize(std::size_t number_of_blocks) const
    {
      return number_of_blocks * row_ids_.size();
    }

    /// @brief Returns the index for a particular element in the compressed data vector
    /// @param number_of_blocks Total number of block sub-matrices in the overall matrix
    /// @param block Index of the block sub-matrix
    /// @param row Index of the row in the block
    /// @param column Index of the column in the block
    /// @return Index of the element in the compressed data vector
    std::size_t VectorIndex(std::size_t number_of_blocks, std::size_t block, std::size_t row, std::size_t column) const
    {
      if (row >= row_start_.size() - 1 || column >= row_start_.size() - 1 || block >= number_of_blocks)
        throw std::system_error(make_error_code(MicmMatrixErrc::ElementOutOfRange));
      auto begin = std::next(row_ids_.begin(), row_start_[row]);
      auto end = std::next(row_ids_.begin(), row_start_[row + 1]);
      auto elem = std::find(begin, end, column);
      if (elem == end)
        throw std::system_error(make_error_code(MicmMatrixErrc::ZeroElementAccess));
      return std::size_t{ (elem - row_ids_.begin()) + block * row_ids_.size() };
    }

    /// @brief Adds a value to the diagonal of each block in the compressed data vector
    /// @param number_of_blocks Total number of block sub-matrices in the overall matrix
    /// @param data Compressed data vector
    /// @param value Value to add to the diagonal
    void AddToDiagonal(const std::size_t number_of_blocks, auto& data, auto value) const
    {
      for (std::size_t block_start = 0; block_start < number_of_blocks * row_ids_.size(); block_start += row_ids_.size())
        for (const auto& i : diagonal_ids_)
          data[block_start + i] += value;
    }

    /// @brief Returns the indices along the diagonal of each block
    /// @param number_of_blocks Number of block sub-matrices in the overall matrix
    /// @param block_id Block index
    /// @return Vector of indices of non-zero diagonal elements
    std::vector<std::size_t> DiagonalIndices(const std::size_t number_of_blocks, const std::size_t block_id) const
    {
      std::vector<std::size_t> indices;
      indices.reserve(row_start_.size() - 1);
      for (std::size_t i = 0; i < row_start_.size() - 1; ++i)
        if (!IsZero(i, i))
          indices.push_back(VectorIndex(number_of_blocks, block_id, i, i));
      return indices;
    }

   public:
    /// @brief A block-local temporary variable with its own storage
    /// For standard ordering: single value
    template<typename T>
    class BlockVariable
    {
     public:
      BlockVariable() = default;
      
      T& Get() { return storage_; }
      const T& Get() const { return storage_; }
      
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
      std::size_t group_;
      std::size_t num_blocks_in_group_;  // Always 1 for standard ordering

      /// @brief Get a const element reference for a specific block in this group
      template<typename Arg>
      [[gnu::always_inline]]
      inline decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg) const
      {
        // For standard ordering, group == block
        std::size_t block = group_;
        
        // Check if Arg has RowIndex() method (ConstBlockView from sparse matrix)
        if constexpr (requires { arg.RowIndex(); arg.ColumnIndex(); })
        {
          // It's a ConstBlockView type from a sparse matrix
          auto* source_matrix = arg.GetMatrix();
          return source_matrix->AsVector()[source_matrix->VectorIndex(block, arg.RowIndex(), arg.ColumnIndex())];
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
            return source_matrix->AsVector()[(row_group * source_matrix->NumColumns() + arg.ColumnIndex()) * L + row_in_group];
          }
          else
          {
            // Standard Matrix layout: data_[row * num_cols + col]
            return source_matrix->AsVector()[block * source_matrix->NumColumns() + arg.ColumnIndex()];
          }
        }
        else if constexpr (requires { arg.Get(); })
        {
          // It's a BlockVariable, for standard ordering always single value
          return arg.Get();
        }
        else
        {
          // Unknown type, just return it
          return arg;
        }
      }

     public:
      ConstGroupView(const SparseMatrixType& matrix, std::size_t group)
          : matrix_(matrix), group_(group), num_blocks_in_group_(1)
      {
      }

      auto GetConstBlockView(std::size_t row, std::size_t col) const
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

      std::size_t NumberOfBlocks() const { return matrix_.NumberOfBlocks(); }
      std::size_t NumRows() const { return matrix_.NumRows(); }
      std::size_t NumColumns() const { return matrix_.NumColumns(); }
    };

    /// @brief GroupView provides a view of a single group of blocks for iteration
    /// For standard ordering: L=1, so each group contains 1 block
    template<typename SparseMatrixType>
    class GroupView
    {
     private:
      SparseMatrixType& matrix_;
      std::size_t group_;
      std::size_t num_blocks_in_group_;  // Always 1 for standard ordering

      /// @brief Get an element reference for a specific block in this group
      template<typename Arg>
      [[gnu::always_inline]]
      inline decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg)
      {
        // For standard ordering, group == block
        std::size_t block = group_;
        
        // Check if Arg has RowIndex() method (BlockView from sparse matrix)
        if constexpr (requires { arg.RowIndex(); arg.ColumnIndex(); })
        {
          // It's a BlockView type from a sparse matrix
          auto* source_matrix = arg.GetMatrix();
          return source_matrix->AsVector()[source_matrix->VectorIndex(block, arg.RowIndex(), arg.ColumnIndex())];
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
            return source_matrix->AsVector()[(row_group * source_matrix->NumColumns() + arg.ColumnIndex()) * L + row_in_group];
          }
          else
          {
            // Standard Matrix layout: data_[row * num_cols + col]
            return source_matrix->AsVector()[block * source_matrix->NumColumns() + arg.ColumnIndex()];
          }
        }
        else if constexpr (requires { arg.Get(); })
        {
          // It's a BlockVariable, for standard ordering always single value
          return arg.Get();
        }
        else
        {
          // Unknown type, just return it
          return arg;
        }
      }

     public:
      GroupView(SparseMatrixType& matrix, std::size_t group)
          : matrix_(matrix), group_(group), num_blocks_in_group_(1)
      {
      }

      auto GetConstBlockView(std::size_t row, std::size_t col) const
      {
        return matrix_.GetConstBlockView(row, col);
      }

      auto GetBlockView(std::size_t row, std::size_t col)
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

      std::size_t NumberOfBlocks() const { return matrix_.NumberOfBlocks(); }
      std::size_t NumRows() const { return matrix_.NumRows(); }
      std::size_t NumColumns() const { return matrix_.NumColumns(); }
    };

    /// @brief Returns the number of blocks included in each group of blocks
    /// @return Number of blocks in each group (1 for standard ordering)
    static constexpr std::size_t GroupVectorSize()
    {
      return 1;
    }

    /// @brief Returns the size of each group of blocks in the compressed data vector
    /// @return Size of each group of blocks
    std::size_t GroupSize(std::size_t number_of_non_zero_elements) const
    {
      return number_of_non_zero_elements;
    }

    /// @brief Returns the total number of groups of blocks in the compressed data
    /// @param number_of_blocks Total number of block sub-matrices in the overall matrix
    /// @return Number of groups of blocks (equal to number_of_blocks for standard ordering)
    std::size_t NumberOfGroups(std::size_t number_of_blocks) const
    {
      return number_of_blocks;
    }

   private:
    /// @brief Returns the row ids of each non-zero element in a block
    /// @param block_size Number of rows or columns for each block
    /// @param non_zero_elements Set of non-zero elements in the matrix
    static std::vector<std::size_t> RowIdsVector(
        const std::size_t block_size,
        const std::set<std::pair<std::size_t, std::size_t>> non_zero_elements)
    {
      std::vector<std::size_t> ids;
      ids.reserve(non_zero_elements.size());
      std::transform(
          non_zero_elements.begin(),
          non_zero_elements.end(),
          std::back_inserter(ids),
          [](const std::pair<std::size_t, std::size_t>& elem) { return elem.second; });
      return ids;
    }

    /// @brief Returns the start and end indices of each row in a block in row_ids_
    /// @param block_size Number of rows or columns for each block
    /// @param non_zero_elements Set of non-zero elements in the matrix
    static std::vector<std::size_t> RowStartVector(
        const std::size_t block_size,
        const std::set<std::pair<std::size_t, std::size_t>> non_zero_elements)
    {
      std::vector<std::size_t> starts(block_size + 1, 0);
      std::size_t total_elem = 0;
      std::size_t curr_row = 0;
      for (auto& elem : non_zero_elements)
      {
        while (curr_row < elem.first)
          starts[(curr_row++) + 1] = total_elem;
        ++total_elem;
      }
      // Fill all remaining entries from curr_row + 1 to block_size
      for (std::size_t i = curr_row + 1; i <= block_size; ++i)
        starts[i] = total_elem;
      return starts;
    }

   public:
    /// @brief Returns whether a particular element is always zero
    /// @param row Row index
    /// @param column Column index
    /// @return true if the element is always zero, false otherwise
    bool IsZero(std::size_t row, std::size_t column) const
    {
      if (row >= row_start_.size() - 1 || column >= row_start_.size() - 1)
        throw std::system_error(make_error_code(MicmMatrixErrc::ElementOutOfRange));
      auto begin = std::next(row_ids_.begin(), row_start_[row]);
      auto end = std::next(row_ids_.begin(), row_start_[row + 1]);
      auto elem = std::find(begin, end, column);
      return (elem == end);
    }
  };
}  // namespace micm
