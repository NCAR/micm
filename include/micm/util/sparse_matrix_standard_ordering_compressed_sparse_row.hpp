// Copyright (C) 2023-2025 National Science Foundation-National Center for Atmospheric Research
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
      starts[curr_row + 1] = total_elem;
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
