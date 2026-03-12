// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/matrix_error.hpp>
#include <micm/util/view_category.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace micm
{

  /// Concept for vectorizable matrices
  template<typename T>
  concept VectorizableDense = requires(T t) {
    t.GroupSize();
    t.GroupVectorSize();
    t.NumberOfGroups();
  };

  /// @brief A 2D array class with contiguous memory
  template<class T = double>
  class Matrix
  {
   public:
    // Diagonal markowitz reordering requires an int argument, make sure one is always accessible
    using IntMatrix = Matrix<int>;
    using value_type = T;

    /// @brief A lightweight descriptor for a const column in a matrix
    class ConstColumnView
    {
      friend class Matrix;
      const Matrix* matrix_;
      std::size_t column_index_;
      
      explicit ConstColumnView(const Matrix* matrix, std::size_t column_index)
          : matrix_(matrix),
            column_index_(column_index)
      {
      }

     public:
      using category = DenseMatrixColumnViewTag;
      std::size_t ColumnIndex() const { return column_index_; }
      const Matrix* GetMatrix() const { return matrix_; }
    };

    /// @brief A lightweight descriptor for a mutable column in a matrix
    class ColumnView
    {
      friend class Matrix;
      Matrix* matrix_;
      std::size_t column_index_;
      
      explicit ColumnView(Matrix* matrix, std::size_t column_index)
          : matrix_(matrix),
            column_index_(column_index)
      {
      }

     public:
      using category = DenseMatrixColumnViewTag;
      std::size_t ColumnIndex() const { return column_index_; }
      Matrix* GetMatrix() { return matrix_; }
    };

    /// @brief A row-local temporary variable with its own storage
    class RowVariable
    {
      friend class Matrix;
      T storage_;  // Stack-allocated single value
      
     public:
      using category = BlockVariableTag;
      RowVariable() = default;
      T& Get() { return storage_; }
      const T& Get() const { return storage_; }
    };

   private:
    std::vector<T> data_;
    std::size_t x_dim_;
    std::size_t y_dim_;

    friend class Proxy;
    friend class ConstProxy;
    template<typename, typename> friend class SparseMatrix;

    class Proxy
    {
      Matrix &matrix_;
      std::size_t offset_;
      std::size_t y_dim_;

     public:
      Proxy(Matrix &matrix, std::size_t offset, std::size_t y_dim)
          : matrix_(matrix),
            offset_(offset),
            y_dim_(y_dim)
      {
      }

      Proxy &operator=(const std::vector<T> &other)
      {
        // check that this row matches the expected rectangular matrix dimensions
        if (other.size() < y_dim_)
        {
          std::string msg = "In matrix row assignment from std::vector. Got " + std::to_string(other.size()) +
                            " elements, but expected " + std::to_string(y_dim_);
          throw std::system_error(make_error_code(MicmMatrixErrc::RowSizeMismatch), msg);
        }
        auto other_elem = other.begin();
        for (auto &elem : *this)
        {
          elem = *(other_elem++);
        }
        return *this;
      }
      operator std::vector<T>() const
      {
        return std::vector<T>(this->begin(), this->end());
      }
      std::size_t Size() const
      {
        return y_dim_;
      }
      typename std::vector<T>::iterator begin() noexcept
      {
        return std::next(matrix_.data_.begin(), offset_);
      }
      typename std::vector<T>::const_iterator begin() const noexcept
      {
        return std::next(matrix_.data_.cbegin(), offset_);
      }
      typename std::vector<T>::iterator end() noexcept
      {
        return std::next(matrix_.data_.begin(), offset_ + y_dim_);
      }
      typename std::vector<T>::const_iterator end() const noexcept
      {
        return std::next(matrix_.data_.begin(), offset_ + y_dim_);
      }
      T &operator[](std::size_t y)
      {
        return matrix_.data_[offset_ + y];
      }
    };

    class ConstProxy
    {
      const Matrix &matrix_;
      std::size_t offset_;
      std::size_t y_dim_;

     public:
      ConstProxy(const Matrix &matrix, std::size_t offset, std::size_t y_dim)
          : matrix_(matrix),
            offset_(offset),
            y_dim_(y_dim)
      {
      }
      operator std::vector<T>() const
      {
        return std::vector<T>(this->begin(), this->end());
      }
      std::size_t Size() const
      {
        return y_dim_;
      }
      typename std::vector<T>::const_iterator begin() const noexcept
      {
        return std::next(matrix_.data_.cbegin(), offset_);
      }
      typename std::vector<T>::const_iterator end() const noexcept
      {
        return std::next(matrix_.data_.begin(), offset_ + y_dim_);
      }
      const T &operator[](std::size_t y) const
      {
        return matrix_.data_[offset_ + y];
      }
    };

   public:
    Matrix()
        : x_dim_(0),
          y_dim_(0),
          data_()
    {
    }

    Matrix(std::size_t x_dim, std::size_t y_dim)
        : x_dim_(x_dim),
          y_dim_(y_dim),
          data_(x_dim * y_dim)
    {
    }

    Matrix(std::size_t x_dim, std::size_t y_dim, T initial_value)
        : x_dim_(x_dim),
          y_dim_(y_dim),
          data_(x_dim * y_dim, initial_value)
    {
    }

    Matrix(const std::vector<std::vector<T>> &other)
        : x_dim_(other.size()),
          y_dim_(other.size() == 0 ? 0 : other[0].size()),
          data_(
              [&]() -> std::vector<T>
              {
                std::size_t x_dim = other.size();
                if (x_dim == 0)
                  return std::vector<T>(0);
                std::size_t y_dim = other[0].size();
                std::vector<T> data(x_dim * y_dim);
                auto elem = data.begin();
                for (std::size_t x{}; x < x_dim; ++x)
                {
                  // check that this row matches the expected rectangular matrix dimensions
                  if (other[x].size() != y_dim)
                  {
                    std::string msg = "In matrix constructor from std::vector<std::vector>. Got " +
                                      std::to_string(other[x].size()) + " columns, but expected " + std::to_string(y_dim);
                    throw std::system_error(make_error_code(MicmMatrixErrc::InvalidVector), "");
                  }
                  for (std::size_t y{}; y < y_dim; ++y)
                  {
                    *(elem++) = other[x][y];
                  }
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
      return y_dim_;
    }

    /// @brief Get the number of elements in the underlying vector between
    ///        adjacent columns for the same row
    /// @return The number of elements in the underlying vector between
    ///         adjacent columns for the same row
    std::size_t ColumnStride() const
    {
      return 1;
    }

    /// @brief Set every matrix element to a given value
    /// @param val Value to set each element to
    void Fill(T val)
    {
      std::fill(data_.begin(), data_.end(), val);
    }

    ConstProxy operator[](std::size_t x) const
    {
      return ConstProxy(*this, x * y_dim_, y_dim_);
    }

    Proxy operator[](std::size_t x)
    {
      return Proxy(*this, x * y_dim_, y_dim_);
    }

    Matrix &operator=(T val)
    {
      std::transform(data_.begin(), data_.end(), data_.begin(), [&](auto &_) { return val; });
      return *this;
    }

    /// @brief For each element in the Matrix x and y, perform y = alpha * x + y,
    ///        where alpha is a scalar constant.
    /// @param alpha The scaling scalar to apply to the Matrix x
    /// @param x The input Matrix
    void Axpy(const double &alpha, const Matrix &x)
    {
      auto x_iter = x.AsVector().begin();
      for (auto &y : data_)
        y += alpha * (*(x_iter++));
    }

    /// @brief For each element of the matrix, perform y = max(y, x), where x is a scalar constant
    /// @param x The scalar constant to compare against
    void Max(const T &x)
    {
      for (auto &y : data_)
        y = std::max(y, x);
    }

    /// @brief For each element of the matrix, perform y = min(y, x), where x is a scalar constant
    /// @param x The scalar constant to compare against
    void Min(const T &x)
    {
      for (auto &y : data_)
        y = std::min(y, x);
    }

    void ForEach(const std::function<void(T &, const T &)> f, const Matrix &a)
    {
      auto a_iter = a.AsVector().begin();
      for (auto &elem : data_)
        f(elem, *(a_iter++));
    }

    void ForEach(const std::function<void(T &, const T &, const T &)> f, const Matrix &a, const Matrix &b)
    {
      auto a_iter = a.AsVector().begin();
      auto b_iter = b.AsVector().begin();
      for (auto &elem : data_)
        f(elem, *(a_iter++), *(b_iter++));
    }

    // Copy the values from the other matrix into this one
    void Copy(const Matrix &other)
    {
      if (other.AsVector().size() != this->data_.size())
        throw std::runtime_error("Both matrices must have the same size.");
      this->data_.assign(other.AsVector().begin(), other.AsVector().end());
    }

    void Swap(Matrix &other)
    {
      if (other.AsVector().size() != this->data_.size())
        throw std::runtime_error("Both matrices must have the same size.");
      data_.swap(other.AsVector());
    }

    // Print the matrix to the output stream
    friend std::ostream &operator<<(std::ostream &os, const Matrix &matrix)
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

    // add begin and end iterators
    typename std::vector<T>::iterator begin() noexcept
    {
      return data_.begin();
    }

    typename std::vector<T>::const_iterator begin() const noexcept
    {
      return data_.cbegin();
    }

    typename std::vector<T>::iterator end() noexcept
    {
      return data_.end();
    }

    typename std::vector<T>::const_iterator end() const noexcept
    {
      return data_.cend();
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
      return RowVariable();
    }

    /// @brief Get a row variable with persistent storage for temporary values (const version)
    /// @return A RowVariable with stack-allocated storage
    RowVariable GetRowVariable() const
    {
      return RowVariable();
    }

    /// @brief Apply a function to each row of the matrix
    /// @tparam Func The lambda/function type
    /// @tparam Args The types of the column view arguments
    /// @param func The function to apply to each row
    /// @param args Column views or row variables
    template<typename Func, typename... Args>
    void ForEachRow(Func&& func, Args&&... args)
    {
      for (std::size_t row = 0; row < x_dim_; ++row)
      {
        func(GetRowElement(row, args)...);
      }
    }

    /// @brief Apply a function to each row of the matrix (const version)
    /// @tparam Func The lambda/function type
    /// @tparam Args The types of the column view arguments
    /// @param func The function to apply to each row
    /// @param args Column views or row variables
    template<typename Func, typename... Args>
    void ForEachRow(Func&& func, Args&&... args) const
    {
      for (std::size_t row = 0; row < x_dim_; ++row)
      {
        func(GetRowElement(row, args)...);
      }
    }

    /// @brief ConstGroupView provides a const view of a single row (group of size 1) for iteration
    class ConstGroupView
    {
     private:
      const Matrix& matrix_;
      std::size_t row_;

      /// @brief Get a const element reference for the current row in this group (ColumnView)
      template<DenseMatrixColumnView Arg>
      [[gnu::always_inline]]
      inline decltype(auto) GetRowElement(Arg&& arg) const
      {
        auto* source_matrix = arg.GetMatrix();
        return source_matrix->data_[row_ * source_matrix->y_dim_ + arg.ColumnIndex()];
      }

      /// @brief Get a const element reference for the current row in this group (RowVariable)
      template<BlockVariableView Arg>
      [[gnu::always_inline]]
      inline decltype(auto) GetRowElement(Arg&& arg) const
      {
        return arg.Get();
      }

      /// @brief Get a const element reference for the current row in this group (Vector-like)
      template<VectorLike Arg>
      [[gnu::always_inline]]
      inline decltype(auto) GetRowElement(Arg&& arg) const
      {
        return arg[row_];
      }

     public:
      ConstGroupView(const Matrix& matrix, std::size_t row)
          : matrix_(matrix), row_(row)
      {
      }

      auto GetConstColumnView(std::size_t column_index) const
      {
        return matrix_.GetConstColumnView(column_index);
      }

      RowVariable GetRowVariable() const
      {
        // Stack-allocated single value
        return RowVariable();
      }

      template<typename Func, typename... Args>
      void ForEachRow(Func&& func, Args&&... args) const
      {
        // For Matrix with L=1, just process the single row (no loop needed)
        func(GetRowElement(std::forward<Args>(args))...);
      }

      std::size_t NumRows() const { return matrix_.NumRows(); }
      std::size_t NumColumns() const { return matrix_.NumColumns(); }
    };

    /// @brief GroupView provides a view of a single row (group of size 1) for iteration
    class GroupView
    {
     private:
      Matrix& matrix_;
      std::size_t row_;

      /// @brief Get an element reference for the current row in this group (ColumnView)
      template<DenseMatrixColumnView Arg>
      [[gnu::always_inline]]
      inline decltype(auto) GetRowElement(Arg&& arg)
      {
        auto* source_matrix = arg.GetMatrix();
        return source_matrix->data_[row_ * source_matrix->y_dim_ + arg.ColumnIndex()];
      }

      /// @brief Get an element reference for the current row in this group (RowVariable)
      template<BlockVariableView Arg>
      [[gnu::always_inline]]
      inline decltype(auto) GetRowElement(Arg&& arg)
      {
        return arg.Get();
      }

      /// @brief Get an element reference for the current row in this group (Vector-like)
      template<VectorLike Arg>
      [[gnu::always_inline]]
      inline decltype(auto) GetRowElement(Arg&& arg)
      {
        return arg[row_];
      }

     public:
      GroupView(Matrix& matrix, std::size_t row)
          : matrix_(matrix), row_(row)
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
        // Stack-allocated single value
        return RowVariable();
      }

      template<typename Func, typename... Args>
      void ForEachRow(Func&& func, Args&&... args)
      {
        // For Matrix with L=1, just process the single row (no loop needed)
        func(GetRowElement(std::forward<Args>(args))...);
      }

      std::size_t NumRows() const { return matrix_.NumRows(); }
      std::size_t NumColumns() const { return matrix_.NumColumns(); }
    };

    /// @brief Create a function that can be applied to matrices and vectors
    /// 
    /// Creates a reusable callable that validates matrix dimensions and applies a user function
    /// row-by-row. For standard Matrix (L=1), each row is processed individually.
    /// 
    /// @tparam Func The lambda/function type
    /// @tparam Args The matrix and vector types
    /// @param func The function to wrap - receives GroupView objects for matrices and vectors
    /// @param args The matrices and vectors to validate and capture dimensions from
    /// @return A callable that validates dimensions and applies the function
    /// 
    /// @note Validation occurs in two phases:
    ///   1. At function creation: Validates row counts match across all matrices and vector sizes
    ///   2. At invocation: Re-validates dimensions in case matrices/vectors were resized
    /// 
    /// @note Column view creation happens inside user lambda and is validated
    ///       at invocation time, not at function creation time. Ensure all column indices
    ///       are within matrix bounds to avoid runtime errors.
    /// 
    /// @throws std::system_error if column counts don't match at creation, vectors have wrong sizes
    ///         at creation, or if at invocation time: matrices/vectors have mismatched row counts,
    ///         column counts don't match creation, or column indices are out of bounds
    template<typename Func, typename... Args>
    static auto Function(Func&& func, Args&... args)
    {
      // Capture column counts for matrices at creation time using helper
      // Row counts can differ between args at creation, but must match at invocation
      auto populate_cols = [](auto&... args_inner) {
        std::vector<std::size_t> cols(sizeof...(args_inner));
        std::size_t idx = 0;
        ([&](auto& arg) {
          using ArgType = std::remove_cvref_t<decltype(arg)>;
          if constexpr (VectorLike<ArgType>) {
            cols[idx] = 0;  // Not used for vectors
          }
          else {
            cols[idx] = arg.NumColumns();
          }
          ++idx;
        }(args_inner), ...);
        return cols;
      };
      
      std::vector<std::size_t> num_cols = populate_cols(args...);

      // Store in variable to ensure fold expression completes before lambda construction
      auto result = [func = std::forward<Func>(func), num_cols = std::move(num_cols)](auto&&... invoked_args) mutable {
        // Validate dimensions and determine row count in a single pass
        std::size_t num_rows = 0;
        bool found_first = false;
        std::size_t idx = 0;
        
        ([&](auto& arg) {
          using ArgType = std::remove_cvref_t<decltype(arg)>;
          
          if constexpr (VectorLike<ArgType>) {
            // Vector - validate size matches row count
            if (!found_first)
            {
              num_rows = arg.size();
              found_first = true;
            }
            else if (arg.size() != num_rows)
            {
              throw std::system_error(
                  make_error_code(MicmMatrixErrc::InvalidVector),
                  "Vector size must match matrix row count. Expected " + std::to_string(num_rows) + 
                      " elements but got " + std::to_string(arg.size()));
            }
          }
          else if constexpr (requires { arg.NumRows(); arg.NumColumns(); }) {
            // Matrix - validate dimensions
            if (!found_first)
            {
              num_rows = arg.NumRows();
              found_first = true;
            }
            else
            {
              if (arg.NumRows() != num_rows)
              {
                throw std::system_error(
                    make_error_code(MicmMatrixErrc::InvalidVector),
                    "All matrices must have the same number of rows when invoking function. Expected " + 
                        std::to_string(num_rows) + " rows but got " + std::to_string(arg.NumRows()));
              }
            }
            
            // Always validate column count against captured value
            if (arg.NumColumns() != num_cols[idx])
            {
              throw std::system_error(
                  make_error_code(MicmMatrixErrc::InvalidVector),
                  "Matrix column count does not match. Expected " + std::to_string(num_cols[idx]) + 
                      " columns but got " + std::to_string(arg.NumColumns()));
            }
          }
          ++idx;
        }(invoked_args), ...);
        
        // Iterate over rows, treating each row as a group of size 1
        for (std::size_t row = 0; row < num_rows; ++row)
        {
          // Use ConstGroupView if matrix is const, otherwise use GroupView
          // For vectors, just pass them through
          func([&](auto&& arg) -> decltype(auto) {
            using ArgType = std::remove_reference_t<decltype(arg)>;
            using ArgTypeNoConst = std::remove_const_t<ArgType>;
            if constexpr (VectorLike<std::remove_cvref_t<ArgType>>)
            {
              // Vector: just forward it
              return std::forward<decltype(arg)>(arg);
            }
            else
            {
              // Matrix: create appropriate GroupView
              if constexpr (std::is_const_v<ArgType>)
              {
                return typename ArgTypeNoConst::ConstGroupView(arg, row);
              }
              else
              {
                return typename ArgTypeNoConst::GroupView(arg, row);
              }
            }
          }(invoked_args)...);
        }
      };
      return result;
    }

   private:
    /// @brief Get an element reference for a row (ColumnView)
    template<DenseMatrixColumnView Arg>
    [[gnu::always_inline]]
    inline decltype(auto) GetRowElement(std::size_t row, Arg&& arg)
    {
      auto* source_matrix = arg.GetMatrix();
      return source_matrix->data_[row * source_matrix->y_dim_ + arg.ColumnIndex()];
    }

    /// @brief Get an element reference for a row (RowVariable)
    template<BlockVariableView Arg>
    [[gnu::always_inline]]
    inline decltype(auto) GetRowElement(std::size_t row, Arg&& arg)
    {
      return arg.Get();
    }

    /// @brief Get an element reference for a row (Vector-like)
    template<VectorLike Arg>
    [[gnu::always_inline]]
    inline decltype(auto) GetRowElement(std::size_t row, Arg&& arg)
    {
      return arg[row];
    }

    /// @brief Get a const element reference for a row (ColumnView) - const version
    template<DenseMatrixColumnView Arg>
    [[gnu::always_inline]]
    inline decltype(auto) GetRowElement(std::size_t row, Arg&& arg) const
    {
      auto* source_matrix = arg.GetMatrix();
      return source_matrix->data_[row * source_matrix->y_dim_ + arg.ColumnIndex()];
    }

    /// @brief Get a const element reference for a row (RowVariable) - const version
    template<BlockVariableView Arg>
    [[gnu::always_inline]]
    inline decltype(auto) GetRowElement(std::size_t row, Arg&& arg) const
    {
      return arg.Get();
    }

    /// @brief Get a const element reference for a row (Vector-like) - const version
    template<VectorLike Arg>
    [[gnu::always_inline]]
    inline decltype(auto) GetRowElement(std::size_t row, Arg&& arg) const
    {
      return arg[row];
    }
  };

  using StandardDenseMatrix = Matrix<double>;

  // ============================================================================
  // Grouping Strategy Specialization
  // ============================================================================

  /// @brief Matrix always uses simple grouping (L==1)
  template<typename T>
  struct GroupingStrategy<Matrix<T>>
  {
    using type = SimpleGroupingTag;
  };

}  // namespace micm
