// Copyright (C) 2025-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/sparse_matrix.hpp>
#include <micm/util/vector_matrix.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace micm
{

  /// @brief Defines the ordering of SparseMatrix object data in Compressed Sparse Column
  ///        format into blocks of rows to encourage vectorization
  ///
  /// Data is stored with sets of blocks in the block diagonal matrix as the highest
  /// level structure, then by column, then by non-zero rows in each column, then by
  /// individual blocks in the set of blocks.
  ///
  /// The template argument is the number of blocks per set of blocks and should be
  /// approximately the size of the vector register.
  template<std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  class SparseMatrixVectorOrderingCompressedSparseColumn
  {
    /// @brief Column ids of each non-zero element in a block
    std::vector<std::size_t> column_ids_;
    /// @brief Start and end indices of each column in a block in column_ids_
    std::vector<std::size_t> column_start_;
    /// @brief Indices along the diagonal of each block
    std::vector<std::size_t> diagonal_ids_;

   protected:
    SparseMatrixVectorOrderingCompressedSparseColumn() = default;

    SparseMatrixVectorOrderingCompressedSparseColumn(
        std::size_t number_of_blocks,
        std::size_t block_size,
        const std::set<std::pair<std::size_t, std::size_t>>& non_zero_elements)
        : column_ids_(ColumnIdsVector(block_size, non_zero_elements)),
          column_start_(ColumnStartVector(block_size, non_zero_elements)),
          diagonal_ids_(DiagonalIndices(number_of_blocks, 0))
    {
    }

    SparseMatrixVectorOrderingCompressedSparseColumn& operator=(
        const std::tuple<std::size_t, std::size_t, std::set<std::pair<std::size_t, std::size_t>>>& number_size_elements)
    {
      column_ids_ = ColumnIdsVector(std::get<1>(number_size_elements), std::get<2>(number_size_elements));
      column_start_ = ColumnStartVector(std::get<1>(number_size_elements), std::get<2>(number_size_elements));
      diagonal_ids_ = DiagonalIndices(std::get<0>(number_size_elements), 0);
      return *this;
    }

    /// @brief Returns the size of the compressed data vector
    /// @param number_of_blocks Number of block sub-matrices in the overall matrix
    /// @return Size of the compressed data vector
    std::size_t VectorSize(std::size_t number_of_blocks) const
    {
      return std::ceil((double)number_of_blocks / (double)L) * L * column_ids_.size();
    }

    /// @brief Returns the index for a particular element in the compressed data vector
    /// @param number_of_blocks Total number of block sub-matrices in the overall matrix
    /// @param block Index of the block sub-matrix
    /// @param row Index of the row in the block
    /// @param column Index of the column in the block
    /// @return Index of the element in the compressed data vector
    std::size_t VectorIndex(std::size_t number_of_blocks, std::size_t block, std::size_t row, std::size_t column) const
    {
      assert(
          column < column_start_.size() - 1 && row < column_start_.size() - 1 && block < number_of_blocks &&
          "element out of range");
      auto begin = std::next(column_ids_.begin(), column_start_[column]);
      auto end = std::next(column_ids_.begin(), column_start_[column + 1]);
      auto elem = std::find(begin, end, row);
      assert(elem != end && "zero element access");
      return std::size_t{ (elem - column_ids_.begin()) * L + block % L + (block / L) * L * column_ids_.size() };
    }

    /// @brief Adds a value to the diagonal of each block in the compressed data vector
    /// @param number_of_blocks Total number of block sub-matrices in the overall matrix
    /// @param data Compressed data vector
    /// @param value Value to add to the diagonal
    void AddToDiagonal(const std::size_t number_of_blocks, auto& data, const auto value) const
    {
      for (std::size_t i_group = 0; i_group < number_of_blocks; i_group += L)
      {
        for (const auto& i : diagonal_ids_)
        {
          auto elem = std::next(data.begin(), i_group * column_ids_.size() + i);
          for (std::size_t i_block = 0; i_block < L; ++i_block)
          {
            elem[i_block] += value;
          }
        }
      }
    }

    /// @brief Convert row and column indices to a block-relative offset
    /// @param row The row index
    /// @param col The column index
    /// @return The block-relative offset (elem_position * L), matching what
    ///         VectorIndex(0, row, col) returns and what GetBlockElement consumes.
    std::size_t VectorIndexFromRowColumn(std::size_t row, std::size_t col) const
    {
      assert(col < column_start_.size() - 1 && row < column_start_.size() - 1 && "element out of range");
      auto begin = std::next(column_ids_.begin(), column_start_[col]);
      auto end = std::next(column_ids_.begin(), column_start_[col + 1]);
      auto elem = std::find(begin, end, row);
      assert(elem != end && "zero element access");
      return std::distance(column_ids_.begin(), elem) * L;
    }

    /// @brief Extract the block-relative offset from a VectorIndex(0, row, col) result
    /// @param vector_index_block_zero The result of VectorIndex(0, row, col)
    /// @return The offset used by GetBlockElement to locate the element within a group.
    ///
    /// For vector ordering this is `elem_position * L` (the raw offset within a
    /// group), NOT `elem_position` (0..number_of_non_zero_elements-1). Returning
    /// the raw offset lets GetBlockElement skip a `* L` per element access, which
    /// balances the `/ L` that a naive implementation would do here. Both cancel
    /// out. The stored value in ConstBlockView / BlockView is opaque to callers,
    /// so this is safe.
    std::size_t ElementPositionFromVectorIndex(std::size_t vector_index_block_zero) const
    {
      return vector_index_block_zero;
    }

    std::vector<std::size_t> DiagonalIndices(const std::size_t number_of_blocks, const std::size_t block_id) const
    {
      std::vector<std::size_t> indices;
      indices.reserve(column_start_.size() - 1);
      for (std::size_t i = 0; i < column_start_.size() - 1; ++i)
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
    /// For vector ordering: array of L values when L>1, single value when L=1
    template<typename T>
    class BlockVariable
    {
     public:
      using category = BlockVariableTag;
      BlockVariable() = default;

      auto& Get()
      {
        return storage_;
      }
      const auto& Get() const
      {
        return storage_;
      }

     private:
      std::conditional_t<(L > 1), std::array<T, L>, T> storage_;
    };

    /// @brief ConstGroupView provides a const view of a single group of blocks for iteration
    /// For vector ordering: each group contains L blocks (except possibly the last group)
    template<typename SparseMatrixType>
    class ConstGroupView
    {
     public:
      using T = typename SparseMatrixType::value_type;

      /// @brief Enriched const block view returned by GetConstBlockView on a ConstGroupView.
      ///        See CSR variant for rationale.
      struct GroupedConstBlockView
      {
        using category = GroupedSparseMatrixBlockViewTag;
        const T* group_base;
        std::size_t block_offset;
      };

     private:
      const SparseMatrixType& matrix_;
      std::size_t group_;
      std::size_t num_blocks_in_group_;  // May be < L for the last group

      /// @brief Get element from sparse matrix ConstBlockView
      template<SparseMatrixBlockView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg) const
      {
        auto* source_matrix = arg.GetMatrix();
        // arg.ElementPosition() is the raw block-relative offset (elem_position * L);
        // see ElementPositionFromVectorIndex above. block_in_group < L is guaranteed
        // by the ForEachBlock loop bound, so (group_ * L + block_in_group) / L == group_
        // and % L == block_in_group; write the address arithmetic without the div/mod.
        std::size_t block_offset = arg.ElementPosition();
        std::size_t num_non_zero = source_matrix->FlatBlockSize();
        return source_matrix->AsVector()[group_ * num_non_zero * L + block_offset + block_in_group];
      }

      /// @brief Get element from GroupedConstBlockView (fast path)
      template<GroupedSparseMatrixBlockView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg) const
      {
        return arg.group_base[arg.block_offset + block_in_group];
      }

      /// @brief Get element from VectorMatrix ConstColumnView (tiered grouping L>1)
      template<DenseMatrixColumnView Arg>
        requires HasTieredGrouping<std::remove_pointer_t<decltype(std::declval<Arg>().GetMatrix())>>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg) const
      {
        auto* source_matrix = arg.GetMatrix();
        // Function() enforces matching L across all matrix args, so we can use L
        // (the sparse ordering's compile-time L) rather than the dense matrix's
        // runtime GroupVectorSize(). block_in_group < L is the ForEachBlock guarantee.
        return source_matrix->AsVector()
            [(group_ * source_matrix->NumColumns() + arg.ColumnIndex()) * L + block_in_group];
      }

      /// @brief Get element from GroupedColumnView (fast path)
      template<GroupedDenseMatrixColumnView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg) const
      {
        return arg.base[block_in_group];
      }

      /// @brief Get element from Matrix or VectorMatrix<1> ConstColumnView (simple grouping L==1)
      template<DenseMatrixColumnView Arg>
        requires HasSimpleGrouping<std::remove_pointer_t<decltype(std::declval<Arg>().GetMatrix())>>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg) const
      {
        // Calculate the actual block index from group and block_in_group
        std::size_t block = group_ * L + block_in_group;
        auto* source_matrix = arg.GetMatrix();

        // Standard Matrix layout: data_[row * num_cols + col]
        return source_matrix->AsVector()[block * source_matrix->NumColumns() + arg.ColumnIndex()];
      }

      /// @brief Get element from BlockVariable (tiered grouping L>1)
      template<BlockVariableView Arg>
        requires(L > 1)
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg) const
      {
        // Vector ordering: BlockVariable has array storage
        return arg.Get()[block_in_group];
      }

      /// @brief Get element from BlockVariable (simple grouping L==1)
      template<BlockVariableView Arg>
        requires(L == 1)
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg) const
      {
        // L=1 case: BlockVariable has single value storage
        return arg.Get();
      }

      /// @brief Get element from Vector-like (tiered grouping L>1)
      template<VectorLike Arg>
        requires(L > 1)
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg) const
      {
        return arg[group_ * L + block_in_group];
      }

      /// @brief Get element from Vector-like (simple grouping L==1)
      template<VectorLike Arg>
        requires(L == 1)
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg) const
      {
        return arg[group_];
      }

     public:
      ConstGroupView(const SparseMatrixType& matrix, std::size_t group)
          : matrix_(matrix),
            group_(group)
      {
        // Calculate how many blocks are in this group
        std::size_t total_blocks = matrix_.NumberOfBlocks();
        std::size_t start_block = group * L;
        num_blocks_in_group_ = std::min(L, total_blocks - start_block);
      }

      /// @brief Returns a grouped const block view whose group base pointer is
      ///        precomputed for this ConstGroupView's group.
      GroupedConstBlockView GetConstBlockView(std::size_t vector_index) const
      {
        return { matrix_.AsVector().data() + group_ * matrix_.FlatBlockSize() * L, vector_index };
      }

      auto GetConstBlockView(std::size_t row, std::size_t col) const
      {
        return matrix_.GetConstBlockView(row, col);
      }

      auto GetBlockVariable() const
      {
        return BlockVariable<T>();
      }

      template<typename Func, typename... Args>
      void ForEachBlock(Func&& func, Args&&... args) const
      {
        // Vector-ordered matrix storage is padded to ceil(N/L)*L cells, so it is always
        // safe to process L blocks per group for matrix args. VectorLike args (e.g.
        // std::vector<T>), however, have exactly N elements and would OOB past the
        // vector's real size, so we fall back to the runtime bound whenever any arg is
        // VectorLike.
        constexpr bool has_vector_arg = (VectorLike<std::remove_cvref_t<Args>> || ...);
        if constexpr (has_vector_arg)
        {
          for (std::size_t block_in_group = 0; block_in_group < num_blocks_in_group_; ++block_in_group)
          {
            func(GetBlockElement(block_in_group, std::forward<Args>(args))...);  // NOLINT(bugprone-use-after-move)
          }
        }
        else
        {
          // Fast path: L is a compile-time constant so the compiler fully unrolls / vectorizes.
          for (std::size_t block_in_group = 0; block_in_group < L; ++block_in_group)
          {
            func(GetBlockElement(block_in_group, std::forward<Args>(args))...);  // NOLINT(bugprone-use-after-move)
          }
        }
      }

      std::size_t NumberOfBlocks() const
      {
        return matrix_.NumberOfBlocks();
      }
      std::size_t NumRows() const
      {
        return matrix_.NumRows();
      }
      std::size_t NumColumns() const
      {
        return matrix_.NumColumns();
      }
    };

    /// @brief GroupView provides a view of a single group of blocks for iteration
    /// For vector ordering: each group contains L blocks (except possibly the last group)
    template<typename SparseMatrixType>
    class GroupView
    {
     public:
      using T = typename SparseMatrixType::value_type;

      /// @brief Enriched mutable block view returned by GetBlockView on a GroupView.
      ///        See ConstGroupView::GroupedConstBlockView for rationale.
      struct GroupedBlockView
      {
        using category = GroupedSparseMatrixBlockViewTag;
        T* group_base;
        std::size_t block_offset;
      };
      /// @brief Const variant, for GetConstBlockView on a mutable GroupView.
      struct GroupedConstBlockView
      {
        using category = GroupedSparseMatrixBlockViewTag;
        const T* group_base;
        std::size_t block_offset;
      };

     private:
      SparseMatrixType& matrix_;
      std::size_t group_;
      std::size_t num_blocks_in_group_;  // May be < L for the last group

      /// @brief Get element from sparse matrix BlockView
      template<SparseMatrixBlockView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg)
      {
        auto* source_matrix = arg.GetMatrix();
        // See ConstGroupView::GetBlockElement for the reasoning.
        std::size_t block_offset = arg.ElementPosition();
        std::size_t num_non_zero = source_matrix->FlatBlockSize();
        return source_matrix->AsVector()[group_ * num_non_zero * L + block_offset + block_in_group];
      }

      /// @brief Get element from GroupedBlockView (fast path)
      template<GroupedSparseMatrixBlockView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg)
      {
        return arg.group_base[arg.block_offset + block_in_group];
      }

      /// @brief Get element from VectorMatrix ColumnView (tiered grouping L>1)
      template<DenseMatrixColumnView Arg>
        requires HasTieredGrouping<std::remove_pointer_t<decltype(std::declval<Arg>().GetMatrix())>>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg)
      {
        auto* source_matrix = arg.GetMatrix();
        // See ConstGroupView::GetBlockElement for the reasoning.
        return source_matrix->AsVector()
            [(group_ * source_matrix->NumColumns() + arg.ColumnIndex()) * L + block_in_group];
      }

      /// @brief Get element from GroupedColumnView (fast path)
      template<GroupedDenseMatrixColumnView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg)
      {
        return arg.base[block_in_group];
      }

      /// @brief Get element from Matrix or VectorMatrix<1> ColumnView (simple grouping L==1)
      template<DenseMatrixColumnView Arg>
        requires HasSimpleGrouping<std::remove_pointer_t<decltype(std::declval<Arg>().GetMatrix())>>
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg)
      {
        // Calculate the actual block index from group and block_in_group
        std::size_t block = group_ * L + block_in_group;
        auto* source_matrix = arg.GetMatrix();

        // Standard Matrix layout: data_[row * num_cols + col]
        return source_matrix->AsVector()[block * source_matrix->NumColumns() + arg.ColumnIndex()];
      }

      /// @brief Get element from BlockVariable (tiered grouping L>1)
      template<BlockVariableView Arg>
        requires(L > 1)
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg)
      {
        // Vector ordering: BlockVariable has array storage
        return arg.Get()[block_in_group];
      }

      /// @brief Get element from BlockVariable (simple grouping L==1)
      template<BlockVariableView Arg>
        requires(L == 1)
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg)
      {
        // L=1 case: BlockVariable has single value storage
        return arg.Get();
      }

      /// @brief Get element from Vector-like (tiered grouping L>1)
      template<VectorLike Arg>
        requires(L > 1)
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg)
      {
        return arg[group_ * L + block_in_group];
      }

      /// @brief Get element from Vector-like (simple grouping L==1)
      template<VectorLike Arg>
        requires(L == 1)
      [[gnu::always_inline]]
      decltype(auto) GetBlockElement(std::size_t block_in_group, Arg&& arg)
      {
        return arg[group_];
      }

     public:
      GroupView(SparseMatrixType& matrix, std::size_t group)
          : matrix_(matrix),
            group_(group)
      {
        // Calculate how many blocks are in this group
        std::size_t total_blocks = matrix_.NumberOfBlocks();
        std::size_t start_block = group * L;
        num_blocks_in_group_ = std::min(L, total_blocks - start_block);
      }

      /// @brief Returns a grouped const block view whose group base pointer is
      ///        precomputed for this GroupView's group.
      GroupedConstBlockView GetConstBlockView(std::size_t vector_index) const
      {
        return { matrix_.AsVector().data() + group_ * matrix_.FlatBlockSize() * L, vector_index };
      }

      auto GetConstBlockView(std::size_t row, std::size_t col) const
      {
        return matrix_.GetConstBlockView(row, col);
      }

      /// @brief Returns a grouped mutable block view whose group base pointer is
      ///        precomputed for this GroupView's group.
      GroupedBlockView GetBlockView(std::size_t vector_index)
      {
        return { matrix_.AsVector().data() + group_ * matrix_.FlatBlockSize() * L, vector_index };
      }

      auto GetBlockView(std::size_t row, std::size_t col)
      {
        return matrix_.GetBlockView(row, col);
      }

      auto GetBlockVariable()
      {
        return BlockVariable<T>();
      }

      template<typename Func, typename... Args>
      void ForEachBlock(Func&& func, Args&&... args)
      {
        // See ConstGroupView::ForEachBlock for rationale.
        constexpr bool has_vector_arg = (VectorLike<std::remove_cvref_t<Args>> || ...);
        if constexpr (has_vector_arg)
        {
          for (std::size_t block_in_group = 0; block_in_group < num_blocks_in_group_; ++block_in_group)
          {
            func(GetBlockElement(block_in_group, std::forward<Args>(args))...);  // NOLINT(bugprone-use-after-move)
          }
        }
        else
        {
          for (std::size_t block_in_group = 0; block_in_group < L; ++block_in_group)
          {
            func(GetBlockElement(block_in_group, std::forward<Args>(args))...);  // NOLINT(bugprone-use-after-move)
          }
        }
      }

      std::size_t NumberOfBlocks() const
      {
        return matrix_.NumberOfBlocks();
      }
      std::size_t NumRows() const
      {
        return matrix_.NumRows();
      }
      std::size_t NumColumns() const
      {
        return matrix_.NumColumns();
      }
    };

    /// @brief Returns the number of blocks included in each group of blocks
    /// @return Number of blocks per group
    static constexpr std::size_t GroupVectorSize()
    {
      return L;
    }

    /// @brief Returns the size of each group of blocks in the compressed data vector
    /// @param number_of_non_zero_elements Number of non-zero elements in the matrix
    /// @return Size of each group of blocks
    std::size_t GroupSize() const
    {
      return L * column_ids_.size();
    }

    /// @brief Returns the total number of groups of blocks in the compressed data
    ///        vector, including any partial groups
    /// @param number_of_blocks Total number of block sub-matrices in the overall matrix
    /// @return Number of groups of blocks
    std::size_t NumberOfGroups(std::size_t number_of_blocks) const
    {
      return std::ceil((double)number_of_blocks / (double)L);
    }

   private:
    /// @brief Returns the column ids of each non-zero element in a block
    /// @param block_size Number of rows or columns in each block
    /// @param non_zero_elements Set of non-zero elements in the matrix
    /// @return Column ids of each non-zero element in a block
    std::vector<std::size_t> ColumnIdsVector(
        const std::size_t block_size,
        const std::set<std::pair<std::size_t, std::size_t>>& non_zero_elements) const
    {
      std::vector<std::size_t> ids;
      ids.reserve(non_zero_elements.size());
      std::set<std::pair<std::size_t, std::size_t>> column_ordered_pairs;
      for (const auto& elem : non_zero_elements)
      {
        column_ordered_pairs.insert(std::make_pair(elem.second, elem.first));
      }
      std::transform(
          column_ordered_pairs.begin(),
          column_ordered_pairs.end(),
          std::back_inserter(ids),
          [](const std::pair<std::size_t, std::size_t>& elem) { return elem.second; });
      return ids;
    }

    /// @brief Returns the start and end indices of each column in a block in column_ids_
    /// @param block_size Number of rows or columns in each block
    /// @param non_zero_elements Set of non-zero elements in the matrix
    /// @return Start and end indices of each column in a block in column_ids_
    std::vector<std::size_t> ColumnStartVector(
        const std::size_t block_size,
        const std::set<std::pair<std::size_t, std::size_t>>& non_zero_elements) const
    {
      std::vector<std::size_t> starts(block_size + 1, 0);
      std::size_t total_elem = 0;
      std::size_t curr_row = 0;
      std::set<std::pair<std::size_t, std::size_t>> column_ordered_pairs;
      for (const auto& elem : non_zero_elements)
      {
        column_ordered_pairs.insert(std::make_pair(elem.second, elem.first));
      }
      for (const auto& elem : column_ordered_pairs)
      {
        while (curr_row < elem.first)
        {
          starts[(curr_row++) + 1] = total_elem;
        }
        ++total_elem;
      }
      // Fill all remaining entries from curr_row + 1 to block_size
      for (std::size_t i = curr_row + 1; i <= block_size; ++i)
      {
        starts[i] = total_elem;
      }
      return starts;
    }

   public:
    /// @brief Returns whether a given row and column index is a zero element
    /// @param row Index of the row
    /// @param column Index of the column
    /// @return True if the element is zero, false otherwise
    bool IsZero(const std::size_t row, const std::size_t column) const
    {
      assert(column < column_start_.size() - 1 && row < column_start_.size() - 1 && "element out of range");
      auto begin = std::next(column_ids_.begin(), column_start_[column]);
      auto end = std::next(column_ids_.begin(), column_start_[column + 1]);
      auto elem = std::find(begin, end, row);
      return (elem == end);
    }
  };

}  // namespace micm