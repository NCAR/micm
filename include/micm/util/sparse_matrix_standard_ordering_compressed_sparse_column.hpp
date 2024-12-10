// Copyright (C) 2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/matrix_error.hpp>

#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <iterator>
#include <utility>
#include <vector>
#include <set>
#include <tuple>

namespace micm
{

  /// @brief Defines the ordering of SparseMatrix object data in Compressed Sparse Column format
  ///
  /// Data is stored with blocks in the block diagonal matrix as the highest
  /// level structure, then by column, then by non-zero rows in each column.
  class SparseMatrixStandardOrderingCompressedSparseColumn
  {
    /// @brief Column ids of each non-zero element in a block
    std::vector<std::size_t> column_ids_;
    /// @brief Start and end indices of each column in a block in column_ids_
    std::vector<std::size_t> column_start_;
    /// @brief Indices along the diagonal of each block
    std::vector<std::size_t> diagonal_ids_;

   protected:
    SparseMatrixStandardOrderingCompressedSparseColumn() = default;

    SparseMatrixStandardOrderingCompressedSparseColumn(
        std::size_t number_of_blocks,
        std::size_t block_size,
        std::set<std::pair<std::size_t, std::size_t>> non_zero_elements)
        : column_ids_(ColumnIdsVector(block_size, non_zero_elements)),
          column_start_(ColumnStartVector(block_size, non_zero_elements)),
          diagonal_ids_(DiagonalIndices(number_of_blocks, 0))
    {
    }

    SparseMatrixStandardOrderingCompressedSparseColumn& operator=(
        const std::tuple<std::size_t, std::size_t, std::set<std::pair<std::size_t, std::size_t>>> number_size_elements)
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
      return number_of_blocks * column_ids_.size();
    }

    /// @brief Returns the index for a particular element in the compressed data vector
    /// @param number_of_blocks Total number of block sub-matrices in the overall matrix
    /// @param block Index of the block sub-matrix
    /// @param row Index of the row in the block
    /// @param column Index of the column in the block
    /// @return Index of the element in the compressed data vector
    std::size_t VectorIndex(
        std::size_t number_of_blocks,
        std::size_t block,
        std::size_t row,
        std::size_t column) const
    {
      if (column >= column_start_.size() - 1 || row >= column_start_.size() - 1 || block >= number_of_blocks)
        throw std::system_error(make_error_code(MicmMatrixErrc::ElementOutOfRange));
      auto begin = std::next(column_ids_.begin(), column_start_[column]);
      auto end = std::next(column_ids_.begin(), column_start_[column + 1]);
      auto elem = std::find(begin, end, row);
      if (elem == end)
        throw std::system_error(make_error_code(MicmMatrixErrc::ZeroElementAccess));
      return std::size_t{ (elem - column_ids_.begin()) + block * column_ids_.size() };
    }

    /// @brief Adds a value to the diagonal of each block in the compressed data vector
    /// @param number_of_blocks Total number of block sub-matrices in the overall matrix
    /// @param data Compressed data vector
    /// @param value Value to add to the diagonal
    void AddToDiagonal(
        const std::size_t number_of_blocks,
        auto& data,
        auto value) const
    {
      for (std::size_t block_start = 0; block_start < number_of_blocks * column_ids_.size(); block_start += column_ids_.size())
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
      indices.reserve(column_start_.size() - 1);
      for (std::size_t i = 0; i < column_start_.size() - 1; ++i)
        if (!IsZero(i, i))
          indices.push_back(VectorIndex(number_of_blocks, block_id, i, i));
      return indices;
    }

   private:

    /// @brief Returns the column ids of each non-zero element in a block
    /// @param block_size Number of rows or columns for each block
    /// @param non_zero_elements Set of non-zero elements in the matrix
    static std::vector<std::size_t> ColumnIdsVector(const std::size_t block_size, const std::set<std::pair<std::size_t, std::size_t>> non_zero_elements)
    {
      std::vector<std::size_t> ids;
      ids.reserve(non_zero_elements.size());
      std::set<std::pair<std::size_t, std::size_t>> column_ordered_pairs;
      for (auto& elem : non_zero_elements)
        column_ordered_pairs.insert(std::make_pair(elem.second, elem.first));
      std::transform(
          column_ordered_pairs.begin(),
          column_ordered_pairs.end(),
          std::back_inserter(ids),
          [](const std::pair<std::size_t, std::size_t>& elem) { return elem.second; });
      return ids;
    }

    /// @brief Returns the start and end indices of each column in a block in column_ids_
    /// @param block_size Number of rows or columns for each block
    /// @param non_zero_elements Set of non-zero elements in the matrix
    static std::vector<std::size_t> ColumnStartVector(const std::size_t block_size, const std::set<std::pair<std::size_t, std::size_t>> non_zero_elements)
    {
      std::vector<std::size_t> starts(block_size + 1, 0);
      std::size_t total_elem = 0;
      std::size_t curr_col = 0;
      std::set<std::pair<std::size_t, std::size_t>> column_ordered_pairs;
      for (auto& elem : non_zero_elements)
        column_ordered_pairs.insert(std::make_pair(elem.second, elem.first));
      for (auto& elem : column_ordered_pairs)
      {
        while (curr_col < elem.first)
          starts[(curr_col++) + 1] = total_elem;
        ++total_elem;
      }
      starts[curr_col + 1] = total_elem;
      return starts;
    }

   public:

    /// @brief Returns whether a particular element is always zero
    /// @param row Row index
    /// @param column Column index
    /// @return true if the element is always zero, false otherwise
    bool IsZero(std::size_t row, std::size_t column) const
    {
      if (column >= column_start_.size() - 1 || row >= column_start_.size() - 1)
        throw std::system_error(make_error_code(MicmMatrixErrc::ElementOutOfRange));
      auto begin = std::next(column_ids_.begin(), column_start_[column]);
      auto end = std::next(column_ids_.begin(), column_start_[column + 1]);
      auto elem = std::find(begin, end, row);
      return (elem == end);
    }
  };
} // namespace micm