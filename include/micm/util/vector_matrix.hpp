// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/matrix_error.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <system_error>
#include <vector>

#ifndef MICM_DEFAULT_VECTOR_SIZE
  #define MICM_DEFAULT_VECTOR_SIZE 4
#endif

namespace micm
{

  /// @brief A 2D array class with contiguous memory structured to encourage vectorization
  ///
  /// The memory layout groups rows into groups whose size can be set such that for a single
  /// column, the group of rows can fit in the vector register.
  ///
  /// The template arguments are the type of the matrix elements and the size of the number
  /// of rows per group.
  template<class T, std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
  class VectorMatrix
  {
   public:
    // Diagonal markowitz reordering requires an int argument, make sure one is always accessible
    using IntMatrix = VectorMatrix<int, L>;
    using value_type = T;

    /// @brief A lightweight descriptor for a const column in a matrix
    class ConstColumnView
    {
      friend class VectorMatrix;
      const VectorMatrix* matrix_;
      std::size_t column_index_;
      
      explicit ConstColumnView(const VectorMatrix* matrix, std::size_t column_index)
          : matrix_(matrix),
            column_index_(column_index)
      {
      }

     public:
      std::size_t ColumnIndex() const { return column_index_; }
      const VectorMatrix* GetMatrix() const { return matrix_; }
    };

    /// @brief A lightweight descriptor for a mutable column in a matrix
    class ColumnView
    {
      friend class VectorMatrix;
      VectorMatrix* matrix_;
      std::size_t column_index_;
      
      explicit ColumnView(VectorMatrix* matrix, std::size_t column_index)
          : matrix_(matrix),
            column_index_(column_index)
      {
      }

     public:
      std::size_t ColumnIndex() const { return column_index_; }
      VectorMatrix* GetMatrix() { return matrix_; }
    };

    /// @brief A row-local temporary variable with its own storage
    class RowVariable
    {
      friend class VectorMatrix;
      alignas(32) std::array<T, L> storage_;  // Stack-allocated SIMD-aligned array
      
     public:
      RowVariable() = default;
      std::array<T, L>& Get() { return storage_; }
      const std::array<T, L>& Get() const { return storage_; }
    };

   private:
   protected:
    alignas(32) std::vector<T> data_;  // SIMD-aligned for vectorization
    std::size_t x_dim_;  // number of rows
    std::size_t y_dim_;  // number of columns

   private:
    friend class Proxy;
    friend class ConstProxy;
    
    // Allow SparseMatrix::GroupView to access data_ for cross-matrix operations
    template<typename U, typename OrderingPolicy>
    friend class SparseMatrix;

    class Proxy
    {
      VectorMatrix &matrix_;
      std::size_t group_index_;
      std::size_t row_index_;
      std::size_t y_dim_;

     public:
      Proxy(VectorMatrix &matrix, std::size_t group_index, std::size_t row_index, std::size_t y_dim)
          : matrix_(matrix),
            group_index_(group_index),
            row_index_(row_index),
            y_dim_(y_dim)
      {
      }

      Proxy &operator=(const std::vector<T> &other)
      {
        if (other.size() < y_dim_)
        {
          std::string msg = "In vector matrix row assignment from std::vector. Got " + std::to_string(other.size()) +
                            " elements, but expected " + std::to_string(y_dim_);
          throw std::system_error(make_error_code(MicmMatrixErrc::RowSizeMismatch), msg);
        }
        auto iter = std::next(matrix_.data_.begin(), group_index_ * y_dim_ * L + row_index_);
        std::for_each(
            other.begin(),
            std::next(other.begin(), y_dim_),
            [&](T const &elem)
            {
              *iter = elem;
              // don't iterate past the end of the vector
              std::size_t remaining_elements = std::distance(iter, matrix_.data_.end());
              iter += std::min(L, remaining_elements);
            });
        return *this;
      }

      operator std::vector<T>() const
      {
        std::vector<T> vec(y_dim_);
        auto iter = std::next(matrix_.data_.begin(), group_index_ * y_dim_ * L + row_index_);
        for (auto &elem : vec)
        {
          elem = *iter;
          // don't iterate past the end of the vector
          std::size_t remaining_elements = std::distance(iter, matrix_.data_.end());
          iter += std::min(L, remaining_elements);
        }
        return vec;
      }

      std::size_t Size() const
      {
        return y_dim_;
      }

      T &operator[](std::size_t y)
      {
        return matrix_.data_[(group_index_ * y_dim_ + y) * L + row_index_];
      }
    };

    class ConstProxy
    {
      const VectorMatrix &matrix_;
      std::size_t group_index_;
      std::size_t row_index_;
      std::size_t y_dim_;

     public:
      ConstProxy(const VectorMatrix &matrix, std::size_t group_index, std::size_t row_index, std::size_t y_dim)
          : matrix_(matrix),
            group_index_(group_index),
            row_index_(row_index),
            y_dim_(y_dim)
      {
      }

      operator std::vector<T>() const
      {
        std::vector<T> vec(y_dim_);
        auto iter = std::next(matrix_.data_.begin(), group_index_ * y_dim_ * L + row_index_);
        for (auto &elem : vec)
        {
          elem = *iter;
          iter += L;
        }
        return vec;
      }

      std::size_t Size() const
      {
        return y_dim_;
      }

      const T &operator[](std::size_t y) const
      {
        return matrix_.data_[(group_index_ * y_dim_ + y) * L + row_index_];
      }
    };

   public:
    VectorMatrix()
        : x_dim_(0),
          y_dim_(0),
          data_()
    {
    }

    VectorMatrix(std::size_t x_dim, std::size_t y_dim)
        : x_dim_(x_dim),
          y_dim_(y_dim),
          data_(std::ceil(x_dim / (double)L) * L * y_dim)
    {
    }

    VectorMatrix(std::size_t x_dim, std::size_t y_dim, T initial_value)
        : x_dim_(x_dim),
          y_dim_(y_dim),
          data_(std::ceil(x_dim / (double)L) * L * y_dim, initial_value)
    {
    }

    VectorMatrix(const std::vector<std::vector<T>> &other)
        : x_dim_(other.size()),
          y_dim_(other.size() == 0 ? 0 : other[0].size()),
          data_(
              [&]() -> std::vector<T>
              {
                std::size_t x_dim = other.size();
                if (x_dim == 0)
                  return std::vector<T>(0);
                std::size_t y_dim = other[0].size();
                std::vector<T> data(std::ceil(x_dim / (double)L) * L * y_dim);
                std::size_t i_row = 0;
                for (auto &other_row : other)
                {
                  if (other_row.size() != y_dim)
                  {
                    std::string msg = "In vector matrix constructor from std::vector<std::vector>. Got " +
                                      std::to_string(other_row.size()) + " columns, but expected " + std::to_string(y_dim);
                    throw std::system_error(make_error_code(MicmMatrixErrc::InvalidVector), msg);
                  }
                  auto iter = std::next(data.begin(), std::floor(i_row / (double)L) * y_dim * L + i_row % L);
                  for (auto &elem : other_row)
                  {
                    *iter = elem;
                    // don't iterate past the end of the vector
                    std::size_t remaining_elements = std::distance(iter, data.end());
                    iter += std::min(L, remaining_elements);
                  }
                  ++i_row;
                }
                return data;
              }())
    {
    }

    std::size_t NumRows() const
    {
      return x_dim_;
    }

    std::size_t NumColumns() const
    {
      return y_dim_;
    }

    /// @brief Get the number of elements in the underlying vector between
    ///        adjacent rows for the same column
    /// @return The number of elements in the underlying vector between
    ///         adjacent rows for the same column
    std::size_t RowStride() const
    {
      return 1;
    }

    /// @brief Get the number of elements in the underlying vector between
    ///        adjacent columns for the same row
    /// @return The number of elements in the underlying vector between
    ///         adjacent columns for the same row
    std::size_t ColumnStride() const
    {
      return L;
    }

    std::size_t NumberOfGroups() const
    {
      return std::ceil(x_dim_ / (double)L);
    }

    std::size_t GroupSize() const
    {
      return L * y_dim_;
    }

    static constexpr std::size_t GroupVectorSize()
    {
      return L;
    }

    void Print() const
    {
      for (std::size_t i = 0; i < x_dim_; ++i)
      {
        for (std::size_t j = 0; j < y_dim_; ++j)
        {
          std::cout << (*this)[i][j] << " ";
        }
        std::cout << std::endl;
      }
    }

    /// @brief Set every matrix element to a given value
    /// @param val Value to set each element to
    void Fill(T val)
    {
      std::fill(data_.begin(), data_.end(), val);
    }

    ConstProxy operator[](std::size_t x) const
    {
      return ConstProxy(*this, std::floor(x / L), x % L, y_dim_);
    }

    Proxy operator[](std::size_t x)
    {
      return Proxy(*this, std::floor(x / L), x % L, y_dim_);
    }

    VectorMatrix &operator=(T val)
    {
      std::transform(data_.begin(), data_.end(), data_.begin(), [&](auto &_) { return val; });
      return *this;
    }

    /// @brief For each element in the VectorMatrix x and y, perform y = alpha * x + y,
    ///        where alpha is a scalar constant.
    /// @param alpha The scaling scalar to apply to the VectorMatrix x
    /// @param x The input VectorMatrix
    void Axpy(const double &alpha, const VectorMatrix &x)
    {
      auto y_iter = data_.begin();
      auto x_iter = x.AsVector().begin();
      const std::size_t n = std::floor(x_dim_ / L) * L * y_dim_;
      for (std::size_t i = 0; i < n; ++i)
      {
        *(y_iter++) += alpha * (*(x_iter++));
      }
      const std::size_t l = x_dim_ % L;
      for (std::size_t i = 0; i < y_dim_; ++i)
      {
        for (std::size_t j = 0; j < l; ++j)
        {
          y_iter[(i * L) + j] += alpha * x_iter[(i * L) + j];
        }
      }
    }

    /// @brief For each element of the VectorMatrix, perform y = max(y, x), where x is a scalar constant
    /// @param x The scalar constant to compare against
    void Max(const T &x)
    {
      for (auto &y : data_)
        y = std::max(y, x);
    }

    /// @brief For each element of the VectorMatrix, perform y = min(y, x), where x is a scalar constant
    /// @param x The scalar constant to compare against
    void Min(const T &x)
    {
      for (auto &y : data_)
        y = std::min(y, x);
    }

    void ForEach(const std::function<void(T &, const T &)> f, const VectorMatrix &a)
    {
      auto this_iter = data_.begin();
      auto a_iter = a.AsVector().begin();
      const std::size_t n = std::floor(x_dim_ / L) * L * y_dim_;
      for (std::size_t i = 0; i < n; ++i)
      {
        f(*(this_iter++), *(a_iter++));
      }
      const std::size_t l = x_dim_ % L;
      for (std::size_t y = 0; y < y_dim_; ++y)
      {
        for (std::size_t x = 0; x < l; ++x)
        {
          f(this_iter[(y * L) + x], a_iter[(y * L) + x]);
        }
      }
    }

    void ForEach(const std::function<void(T &, const T &, const T &)> f, const VectorMatrix &a, const VectorMatrix &b)
    {
      auto this_iter = data_.begin();
      auto a_iter = a.AsVector().begin();
      auto b_iter = b.AsVector().begin();
      const std::size_t n = std::floor(x_dim_ / L) * L * y_dim_;
      for (std::size_t i = 0; i < n; ++i)
      {
        f(*(this_iter++), *(a_iter++), *(b_iter++));
      }
      const std::size_t l = x_dim_ % L;
      if (l > 0)
      {
        for (std::size_t y = 0; y < y_dim_; ++y)
        {
          for (std::size_t x = 0; x < l; ++x)
          {
            f(this_iter[(y * L) + x], a_iter[(y * L) + x], b_iter[(y * L) + x]);
          }
        }
      }
    }

    // Copy the values from the other VectorMatrix into this one
    void Copy(const VectorMatrix &other)
    {
      if (other.AsVector().size() != this->data_.size())
        throw std::runtime_error("Both vector matrices must have the same size.");
      this->data_.assign(other.AsVector().begin(), other.AsVector().end());
    }

    void Swap(VectorMatrix &other)
    {
      if (other.AsVector().size() != this->data_.size())
        throw std::runtime_error("Both vector matrices must have the same size.");
      data_.swap(other.AsVector());
    }

    // Print the VectorMatrix to the output stream
    friend std::ostream &operator<<(std::ostream &os, const VectorMatrix &matrix)
    {
      for (std::size_t i = 0; i < matrix.x_dim_; ++i)
      {
        for (std::size_t j = 0; j < matrix.y_dim_ - 1; ++j)
        {
          os << matrix[i][j] << ',';
        }
        os << matrix[i][matrix.y_dim_ - 1] << std::endl;
      }
      return os;
    }

    std::vector<T> &AsVector()
    {
      return data_;
    }

    const std::vector<T> &AsVector() const
    {
      return data_;
    }

    /// @brief Create a const column view for accessing a column
    /// @param column_index The index of the column
    /// @return A ConstColumnView descriptor
    ConstColumnView GetConstColumnView(std::size_t column_index) const
    {
      if (column_index >= y_dim_)
      {
        throw std::system_error(
            make_error_code(MicmMatrixErrc::ElementOutOfRange),
            "Column index " + std::to_string(column_index) + " out of range for matrix with " +
                std::to_string(y_dim_) + " columns");
      }
      return ConstColumnView(this, column_index);
    }

    /// @brief Create a mutable column view for accessing a column
    /// @param column_index The index of the column
    /// @return A ColumnView descriptor
    ColumnView GetColumnView(std::size_t column_index)
    {
      if (column_index >= y_dim_)
      {
        throw std::system_error(
            make_error_code(MicmMatrixErrc::ElementOutOfRange),
            "Column index " + std::to_string(column_index) + " out of range for matrix with " +
                std::to_string(y_dim_) + " columns");
      }
      return ColumnView(this, column_index);
    }

    /// @brief Get a row variable with persistent storage for temporary values
    /// @return A RowVariable with stack-allocated storage
    RowVariable GetRowVariable()
    {
      // Stack-allocated array of L elements
      return RowVariable();
    }

    /// @brief Apply a function to each row of the matrix (processes L rows at a time)
    /// @tparam Func The lambda/function type
    /// @tparam Args The types of the column view arguments
    /// @param func The function to apply to each row
    /// @param args Column views or row variables
    template<typename Func, typename... Args>
    void ForEachRow(Func&& func, Args&&... args)
    {
      // Process complete groups of L rows
      std::size_t num_groups = std::floor(x_dim_ / (double)L);
      for (std::size_t group = 0; group < num_groups; ++group)
      {
        for (std::size_t row_in_group = 0; row_in_group < L; ++row_in_group)
        {
          std::size_t row = group * L + row_in_group;
          func(GetRowElement(row, group, row_in_group, args)...);
        }
      }
      
      // Process remaining rows (if x_dim_ is not a multiple of L)
      std::size_t remaining = x_dim_ % L;
      if (remaining > 0)
      {
        std::size_t last_group = num_groups;
        for (std::size_t row_in_group = 0; row_in_group < remaining; ++row_in_group)
        {
          std::size_t row = last_group * L + row_in_group;
          func(GetRowElement(row, last_group, row_in_group, args)...);
        }
      }
    }

    /// @brief GroupView provides a view of a single group of L rows for iteration
    class GroupView
    {
     private:
      VectorMatrix& matrix_;
      std::size_t group_;
      std::size_t num_rows_in_group_;  // May be < L for the last group

      /// @brief Get an element reference for a specific row in this group
      template<typename Arg>
      [[gnu::always_inline]]
      inline decltype(auto) GetRowElement(std::size_t row_in_group, Arg&& arg)
      {
        // Check if Arg has GetMatrix() method (ColumnView)
        if constexpr (requires { arg.GetMatrix(); })
        {
          // It's a ColumnView type, access the source matrix's data
          auto* source_matrix = arg.GetMatrix();
          // VectorMatrix layout: data_[(group * y_dim_ + column) * L + row_in_group]
          return source_matrix->data_[(group_ * source_matrix->y_dim_ + arg.ColumnIndex()) * L + row_in_group];
        }
        else if constexpr (requires { arg.Get(); })
        {
          // It's a RowVariable from this GroupView, access the array element
          return arg.Get()[row_in_group];
        }
        else
        {
          // Unknown type, just return it
          return arg;
        }
      }

     public:
      /// @brief Constructor that calculates num_rows_in_group from matrix dimensions
      GroupView(VectorMatrix& matrix, std::size_t group)
          : matrix_(matrix), group_(group)
      {
        // Calculate how many rows are in this group (typically L, except possibly the last group)
        std::size_t total_groups = matrix.NumberOfGroups();
        if (group == total_groups - 1)
        {
          // Last group may have fewer than L rows
          num_rows_in_group_ = matrix.x_dim_ - (total_groups - 1) * L;
        }
        else
        {
          // All other groups have exactly L rows
          num_rows_in_group_ = L;
        }
      }

      /// @brief Constructor with explicit num_rows_in_group
      GroupView(VectorMatrix& matrix, std::size_t group, std::size_t num_rows_in_group)
          : matrix_(matrix), group_(group), num_rows_in_group_(num_rows_in_group)
      {
      }

      auto GetConstColumnView(std::size_t column_index) const
      {
        return matrix_.GetConstColumnView(column_index);
      }

      auto GetColumnView(std::size_t column_index)
      {
        return matrix_.GetColumnView(column_index);
      }

      RowVariable GetRowVariable()
      {
        // Stack-allocated array of L elements
        return RowVariable();
      }

      template<typename Func, typename... Args>
      void ForEachRow(Func&& func, Args&&... args)
      {
        // Tight loop over L rows in this group for vectorization
        for (std::size_t row_in_group = 0; row_in_group < num_rows_in_group_; ++row_in_group)
        {
          func(GetRowElement(row_in_group, std::forward<Args>(args))...);
        }
      }

      std::size_t NumRows() const { return matrix_.NumRows(); }
      std::size_t NumColumns() const { return matrix_.NumColumns(); }
    };

    /// @brief Create a function that can be applied to matrices
    /// @tparam Func The lambda/function type
    /// @tparam Matrices The matrix types
    /// @param func The function to wrap
    /// @param matrices The matrices to validate and capture dimensions from
    /// @return A callable that validates dimensions and applies the function
    template<typename Func, typename... Matrices>
    static auto Function(Func&& func, Matrices&... matrices)
    {
      // Validate that all matrices have the same number of rows
      std::size_t num_rows = 0;
      std::array<std::size_t, sizeof...(Matrices)> num_cols_per_matrix{};
      std::size_t index = 0;
      
      ([&](auto& matrix) {
        if (index == 0)
        {
          num_rows = matrix.NumRows();
        }
        else if (matrix.NumRows() != num_rows)
        {
          throw std::system_error(
              make_error_code(MicmMatrixErrc::InvalidVector),
              "All matrices must have the same number of rows. Expected " + std::to_string(num_rows) +
                  " but got " + std::to_string(matrix.NumRows()));
        }
        num_cols_per_matrix[index] = matrix.NumColumns();
        ++index;
      }(matrices), ...);

      // Return a callable that validates dimensions on invocation and applies the function
      return [func = std::forward<Func>(func), num_rows, num_cols_per_matrix](Matrices&... invoked_matrices) {
        std::size_t idx = 0;
        ([&](auto& matrix) {
          if (matrix.NumRows() != num_rows)
          {
            throw std::system_error(
                make_error_code(MicmMatrixErrc::InvalidVector),
                "Matrix dimensions do not match. Expected " + std::to_string(num_rows) + " rows but got " +
                    std::to_string(matrix.NumRows()));
          }
          if (matrix.NumColumns() != num_cols_per_matrix[idx])
          {
            throw std::system_error(
                make_error_code(MicmMatrixErrc::InvalidVector),
                "Matrix dimensions do not match. Expected " + std::to_string(num_cols_per_matrix[idx]) + 
                    " columns but got " + std::to_string(matrix.NumColumns()));
          }
          ++idx;
        }(invoked_matrices), ...);
        
        // Iterate over groups, processing L rows at a time
        std::size_t num_complete_groups = std::floor(num_rows / (double)L);
        for (std::size_t group = 0; group < num_complete_groups; ++group)
        {
          func(typename std::decay_t<Matrices>::GroupView(invoked_matrices, group, L)...);
        }
        
        // Process remaining rows (if num_rows is not a multiple of L)
        std::size_t remaining = num_rows % L;
        if (remaining > 0)
        {
          func(typename std::decay_t<Matrices>::GroupView(invoked_matrices, num_complete_groups, remaining)...);
        }
      };
    }

   private:
    /// @brief Get an element reference for a row, handling ColumnViews and RowVariables
    template<typename Arg>
    [[gnu::always_inline]]
    inline decltype(auto) GetRowElement(std::size_t row, std::size_t group, std::size_t row_in_group, Arg&& arg)
    {
      // Check if Arg has GetMatrix() method (ColumnView from potentially different matrix)
      if constexpr (requires { arg.GetMatrix(); })
      {
        // It's a ColumnView type, access the source matrix's data
        auto* source_matrix = arg.GetMatrix();
        // VectorMatrix layout: data_[(group * y_dim_ + column) * L + row_in_group]
        return source_matrix->data_[(group * source_matrix->y_dim_ + arg.ColumnIndex()) * L + row_in_group];
      }
      else if constexpr (requires { arg.Get(); })
      {
        // It's a RowVariable, return reference to the array element
        return arg.Get()[row_in_group];
      }
      else
      {
        // Unknown type, just return it
        return arg;
      }
    }
  };

}  // namespace micm
