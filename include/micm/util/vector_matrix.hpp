// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/micm_exception.hpp>
#include <micm/util/padded_vector.hpp>
#include <micm/util/reducers.hpp>
#include <micm/util/scalar_view.hpp>
#include <micm/util/types.hpp>
#include <micm/util/view_category.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

namespace micm
{

  /// @brief A 2D array class with contiguous memory structured to encourage vectorization
  ///
  /// The memory layout groups rows into groups whose size can be set such that for a single
  /// column, the group of rows can fit in the vector register.
  ///
  /// The template arguments are the type of the matrix elements and the size of the number
  /// of rows per group.
  template<class T, Index L = MICM_DEFAULT_VECTOR_SIZE>
  class VectorMatrix
  {
   public:
    // Diagonal markowitz reordering requires an int argument, make sure one is always accessible
    using IntMatrix = VectorMatrix<int, L>;
    using value_type = T;
    class GroupView;
    class ConstGroupView;
    using ViewType = GroupView;
    using ConstViewType = ConstGroupView;
    using HostGroupView = GroupView;
    using ConstHostGroupView = ConstGroupView;
    template<class VecT>
    using VectorType = PaddedVector<VecT, L>;
    template<class ScaT>
    using ScalarType = ScalarView<ScaT>;
    template<class U>
    using SumType = micm::Sum<U>;
    template<class U>
    using MaxType = micm::Max<U>;
    using LOrType = micm::LOr;
    using LAndType = micm::LAnd;

    /// @brief A lightweight descriptor for a const column in a matrix
    class ConstColumnView
    {
      friend class VectorMatrix;
      const VectorMatrix* matrix_;
      Index column_index_;

      explicit ConstColumnView(const VectorMatrix* matrix, Index column_index)
          : matrix_(matrix),
            column_index_(column_index)
      {
      }

     public:
      using category = DenseMatrixColumnViewTag;
      Index ColumnIndex() const
      {
        return column_index_;
      }
      const VectorMatrix* GetMatrix() const
      {
        return matrix_;
      }
    };

    /// @brief A lightweight descriptor for a mutable column in a matrix
    class ColumnView
    {
      friend class VectorMatrix;
      VectorMatrix* matrix_;
      Index column_index_;

      explicit ColumnView(VectorMatrix* matrix, Index column_index)
          : matrix_(matrix),
            column_index_(column_index)
      {
      }

     public:
      using category = DenseMatrixColumnViewTag;
      Index ColumnIndex() const
      {
        return column_index_;
      }
      VectorMatrix* GetMatrix()
      {
        return matrix_;
      }
    };

    /// @brief A row-local temporary variable with its own storage
    class RowVariable
    {
      friend class VectorMatrix;
      alignas(32) std::array<T, L> storage_;  // Stack-allocated SIMD-aligned array

     public:
      using category = BlockVariableTag;
      RowVariable() = default;
      std::array<T, L>& Get()
      {
        return storage_;
      }
      const std::array<T, L>& Get() const
      {
        return storage_;
      }
    };

   protected:
    std::vector<T> data_;  // Memory alignment depends on std::vector's allocator
    Index x_dim_;          // number of rows
    Index y_dim_;          // number of columns

   private:
    friend class Proxy;
    friend class ConstProxy;

    // Allow SparseMatrix::GroupView to access data_ for cross-matrix operations
    template<typename U, typename OrderingPolicy>
    friend class SparseMatrix;

    class Proxy
    {
      VectorMatrix& matrix_;
      Index group_index_;
      Index row_index_;
      Index y_dim_;

     public:
      Proxy(VectorMatrix& matrix, Index group_index, Index row_index, Index y_dim)
          : matrix_(matrix),
            group_index_(group_index),
            row_index_(row_index),
            y_dim_(y_dim)
      {
      }

      Proxy& operator=(const std::vector<T>& other)
      {
        if (other.size() < y_dim_)
        {
          std::string msg = "In vector matrix row assignment from std::vector. Got " + std::to_string(other.size()) +
                            " elements, but expected " + std::to_string(y_dim_);
          throw MicmException(MICM_ERROR_CATEGORY_MATRIX, MICM_MATRIX_ERROR_CODE_ROW_SIZE_MISMATCH, msg);
        }
        auto iter = std::next(matrix_.data_.begin(), group_index_ * y_dim_ * L + row_index_);
        std::for_each(
            other.begin(),
            std::next(other.begin(), y_dim_),
            [&](T const& elem)
            {
              *iter = elem;
              // don't iterate past the end of the vector
              Index remaining_elements = std::distance(iter, matrix_.data_.end());
              iter += std::min(L, remaining_elements);
            });
        return *this;
      }

      operator std::vector<T>() const
      {
        std::vector<T> vec(y_dim_);
        auto iter = std::next(matrix_.data_.begin(), group_index_ * y_dim_ * L + row_index_);
        for (auto& elem : vec)
        {
          elem = *iter;
          // don't iterate past the end of the vector
          Index remaining_elements = std::distance(iter, matrix_.data_.end());
          iter += std::min(L, remaining_elements);
        }
        return vec;
      }

      Index Size() const
      {
        return y_dim_;
      }

      T& operator[](Index y)
      {
        return matrix_.data_[(group_index_ * y_dim_ + y) * L + row_index_];
      }
    };

    class ConstProxy
    {
      const VectorMatrix& matrix_;
      Index group_index_;
      Index row_index_;
      Index y_dim_;

     public:
      ConstProxy(const VectorMatrix& matrix, Index group_index, Index row_index, Index y_dim)
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
        for (auto& elem : vec)
        {
          elem = *iter;
          iter += L;
        }
        return vec;
      }

      Index Size() const
      {
        return y_dim_;
      }

      const T& operator[](Index y) const
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

    VectorMatrix(Index x_dim, Index y_dim)
        : x_dim_(x_dim),
          y_dim_(y_dim),
          data_(std::ceil(x_dim / (double)L) * L * y_dim)
    {
    }

    VectorMatrix(Index x_dim, Index y_dim, T initial_value)
        : x_dim_(x_dim),
          y_dim_(y_dim),
          data_(std::ceil(x_dim / (double)L) * L * y_dim, initial_value)
    {
    }

    VectorMatrix(const std::vector<std::vector<T>>& other)
        : x_dim_(other.size()),
          y_dim_(other.empty() ? 0 : other[0].size()),
          data_(
              [&]() -> std::vector<T>
              {
                Index x_dim = other.size();
                if (x_dim == 0)
                {
                  return std::vector<T>(0);
                }
                Index y_dim = other[0].size();
                std::vector<T> data(std::ceil(x_dim / (double)L) * L * y_dim);
                Index i_row = 0;
                for (auto& other_row : other)
                {
                  if (other_row.size() != y_dim)
                  {
                    std::string msg = "In vector matrix constructor from std::vector<std::vector>. Got " +
                                      std::to_string(other_row.size()) + " columns, but expected " + std::to_string(y_dim);
                    throw MicmException(MICM_ERROR_CATEGORY_MATRIX, MICM_MATRIX_ERROR_CODE_INVALID_VECTOR, msg);
                  }
                  auto iter = std::next(data.begin(), (i_row / L) * y_dim * L + i_row % L);
                  for (auto& elem : other_row)
                  {
                    *iter = elem;
                    // don't iterate past the end of the vector
                    Index remaining_elements = std::distance(iter, data.end());
                    iter += std::min(L, remaining_elements);
                  }
                  ++i_row;
                }
                return data;
              }())
    {
    }

    virtual ~VectorMatrix() = default;

    Index NumRows() const
    {
      return x_dim_;
    }

    Index NumColumns() const
    {
      return y_dim_;
    }

    /// @brief Get the number of elements in the underlying vector between
    ///        adjacent rows for the same column
    /// @return The number of elements in the underlying vector between
    ///         adjacent rows for the same column
    Index RowStride() const
    {
      return 1;
    }

    /// @brief Get the number of elements in the underlying vector between
    ///        adjacent columns for the same row
    /// @return The number of elements in the underlying vector between
    ///         adjacent columns for the same row
    Index ColumnStride() const
    {
      return L;
    }

    Index NumberOfGroups() const
    {
      return std::ceil(x_dim_ / (double)L);
    }

    Index GroupSize() const
    {
      return L * y_dim_;
    }

    static constexpr Index GroupVectorSize()
    {
      return L;
    }

    void Print() const
    {
      for (Index i = 0; i < x_dim_; ++i)
      {
        for (Index j = 0; j < y_dim_; ++j)
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

    /// @brief No-op host-to-device sync hook.
    ///
    /// GPU-backed matrix policies (e.g. KokkosDenseMatrix, CudaDenseMatrix)
    /// override this to copy host data to a device mirror. Defined here as a
    /// no-op so shared MatrixPolicy tests and solver code can call it
    /// unconditionally regardless of which matrix policy is in use.
    void CopyToDevice() const
    {
    }

    /// @brief No-op device-to-host sync hook. See CopyToDevice().
    void CopyToHost() const
    {
    }

    /// @brief Creates a vector usable with this matrix type in Function() lambdas
    /// @param n vector size (excluding padding)
    /// @param init initial value for vector elements
    /// @return vector usable in Function() lambdas
    template<class VecT>
    VectorType<VecT> CompatibleVector(Index n, VecT init = VecT{}) const
    {
      return VectorType<VecT>(n, init);
    }

    /// @brief Creates a scalar usable with this matrix type in Function lambda captures
    /// @param init initial value for scalar
    /// @return scalar usable in Function() lambda captures
    template<class ScaT>
    ScalarType<ScaT> CompatibleScalar(ScaT init = ScaT{}) const
    {
      return ScalarType<ScaT>(init);
    }

    ConstProxy operator[](Index x) const
    {
      return ConstProxy(*this, (x / L), x % L, y_dim_);
    }

    Proxy operator[](Index x)
    {
      return Proxy(*this, (x / L), x % L, y_dim_);
    }

    VectorMatrix& operator=(T val)
    {
      std::transform(data_.begin(), data_.end(), data_.begin(), [&](auto& _) { return val; });
      return *this;
    }

    /// @brief For each element in the VectorMatrix x and y, perform y = alpha * x + y,
    ///        where alpha is a scalar constant.
    /// @param alpha The scaling scalar to apply to the VectorMatrix x
    /// @param x The input VectorMatrix
    void Axpy(const Real& alpha, const VectorMatrix& x)
    {
      auto y_iter = data_.begin();
      auto x_iter = x.AsVector().begin();
      const Index n = (x_dim_ / L) * L * y_dim_;
      for (Index i = 0; i < n; ++i)
      {
        *(y_iter++) += alpha * (*(x_iter++));
      }
      const Index l = x_dim_ % L;
      for (Index i = 0; i < y_dim_; ++i)
      {
        for (Index j = 0; j < l; ++j)
        {
          y_iter[(i * L) + j] += alpha * x_iter[(i * L) + j];
        }
      }
    }

    /// @brief For each element of the VectorMatrix, perform y = max(y, x), where x is a scalar constant
    /// @param x The scalar constant to compare against
    void Max(const T& x)
    {
      for (auto& y : data_)
      {
        y = std::max(y, x);
      }
    }

    /// @brief For each element of the VectorMatrix, perform y = min(y, x), where x is a scalar constant
    /// @param x The scalar constant to compare against
    void Min(const T& x)
    {
      for (auto& y : data_)
      {
        y = std::min(y, x);
      }
    }

    void ForEach(const std::function<void(T&, const T&)>& f, const VectorMatrix& a)
    {
      auto this_iter = data_.begin();
      auto a_iter = a.AsVector().begin();
      const Index n = (x_dim_ / L) * L * y_dim_;
      for (Index i = 0; i < n; ++i)
      {
        f(*(this_iter++), *(a_iter++));
      }
      const Index l = x_dim_ % L;
      for (Index y = 0; y < y_dim_; ++y)
      {
        for (Index x = 0; x < l; ++x)
        {
          f(this_iter[(y * L) + x], a_iter[(y * L) + x]);
        }
      }
    }

    void ForEach(const std::function<void(T&, const T&, const T&)>& f, const VectorMatrix& a, const VectorMatrix& b)
    {
      auto this_iter = data_.begin();
      auto a_iter = a.AsVector().begin();
      auto b_iter = b.AsVector().begin();
      const Index n = (x_dim_ / L) * L * y_dim_;
      for (Index i = 0; i < n; ++i)
      {
        f(*(this_iter++), *(a_iter++), *(b_iter++));
      }
      const Index l = x_dim_ % L;
      if (l > 0)
      {
        for (Index y = 0; y < y_dim_; ++y)
        {
          for (Index x = 0; x < l; ++x)
          {
            f(this_iter[(y * L) + x], a_iter[(y * L) + x], b_iter[(y * L) + x]);
          }
        }
      }
    }

    // Copy the values from the other VectorMatrix into this one
    void Copy(const VectorMatrix& other)
    {
      if (other.AsVector().size() != this->data_.size())
      {
        throw std::runtime_error("Both vector matrices must have the same size.");
      }
      this->data_.assign(other.AsVector().begin(), other.AsVector().end());
    }

    void Swap(VectorMatrix& other)
    {
      if (other.AsVector().size() != this->data_.size())
      {
        throw std::runtime_error("Both vector matrices must have the same size.");
      }
      data_.swap(other.AsVector());
    }

    // Print the VectorMatrix to the output stream
    friend std::ostream& operator<<(std::ostream& os, const VectorMatrix& matrix)
    {
      for (Index i = 0; i < matrix.x_dim_; ++i)
      {
        for (Index j = 0; j < matrix.y_dim_ - 1; ++j)
        {
          os << matrix[i][j] << ',';
        }
        os << matrix[i][matrix.y_dim_ - 1] << std::endl;
      }
      return os;
    }

    std::vector<T>& AsVector()
    {
      return data_;
    }

    const std::vector<T>& AsVector() const
    {
      return data_;
    }

    /// @brief Create a const column view for accessing a column
    /// @param column_index The index of the column
    /// @return A ConstColumnView descriptor
    ConstColumnView GetConstColumnView(Index column_index) const
    {
      assert(column_index < y_dim_ && "column index out of range");
      return ConstColumnView(this, column_index);
    }

    /// @brief Create a mutable column view for accessing a column
    /// @param column_index The index of the column
    /// @return A ColumnView descriptor
    ColumnView GetColumnView(Index column_index) const
    {
      assert(column_index < y_dim_ && "column index out of range");
      return ColumnView(this, column_index);
    }

    /// @brief Get a row variable with persistent storage for temporary values (const version)
    /// @return A RowVariable with stack-allocated storage
    RowVariable GetRowVariable() const
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
      Index num_groups = x_dim_ / L;
      for (Index group = 0; group < num_groups; ++group)
      {
        for (Index row_in_group = 0; row_in_group < L; ++row_in_group)
        {
          Index row = group * L + row_in_group;
          func(GetRowElement(row, group, row_in_group, args)...);
        }
      }

      // Process remaining rows (if x_dim_ is not a multiple of L)
      Index remaining = x_dim_ % L;
      if (remaining > 0)
      {
        Index last_group = num_groups;
        for (Index row_in_group = 0; row_in_group < remaining; ++row_in_group)
        {
          Index row = last_group * L + row_in_group;
          func(GetRowElement(row, last_group, row_in_group, args)...);
        }
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
      // Process complete groups of L rows
      Index num_groups = x_dim_ / L;
      for (Index group = 0; group < num_groups; ++group)
      {
        for (Index row_in_group = 0; row_in_group < L; ++row_in_group)
        {
          Index row = group * L + row_in_group;
          func(GetRowElement(row, group, row_in_group, args)...);
        }
      }

      // Process remaining rows (if x_dim_ is not a multiple of L)
      Index remaining = x_dim_ % L;
      if (remaining > 0)
      {
        Index last_group = num_groups;
        for (Index row_in_group = 0; row_in_group < remaining; ++row_in_group)
        {
          Index row = last_group * L + row_in_group;
          func(GetRowElement(row, last_group, row_in_group, args)...);
        }
      }
    }

    /// @brief ConstGroupView provides a const view of a single group of L rows for iteration
    class ConstGroupView
    {
     public:
      /// @brief Enriched column view returned by GetConstColumnView on a ConstGroupView.
      ///
      /// Carries a precomputed base_ pointer into this ConstGroupView's slice of the
      /// underlying storage. For VectorMatrix, `base_` points at the first row of the
      /// group's L-row block for `column_index`, so element access is
      /// `arg.base_[row_in_group]` (contiguous) instead of recomputing
      /// `(group * y_dim + column) * L + row_in_group` per element.
      struct GroupedConstColumnView
      {
        using category = GroupedDenseMatrixColumnViewTag;
        const T* base_;
      };

     private:
      const VectorMatrix& matrix_;
      Index group_;
      Index num_rows_in_group_;  // May be < L for the last group

      /// @brief Get a const element reference for a specific row in this group (ColumnView)
      template<DenseMatrixColumnView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        auto* source_matrix = arg.GetMatrix();
        // VectorMatrix layout: data_[(group * y_dim_ + column) * L + row_in_group]
        return source_matrix->data_[(group_ * source_matrix->y_dim_ + arg.ColumnIndex()) * L + row_in_group];
      }

      /// @brief Get a const element reference for a specific row in this group (GroupedColumnView)
      /// Fast path: `base_` already points at row 0 of this group's L-row block.
      template<GroupedDenseMatrixColumnView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg.base_[row_in_group];
      }

      /// @brief Get a const element reference for a specific row in this group (BlockVariable).
      ///        BlockVariable::Get() flavor depends on the source policy:
      ///          - Dense RowVariable is `std::array<T, L>` for all L (including L=1).
      ///          - Sparse L>1 BlockVariable is `std::array<T, L>`.
      ///          - Sparse L=1 BlockVariable is a scalar T.
      ///        Dispatch on whether the returned storage is subscriptable rather
      ///        than on L, so both flavors work for either L.
      template<BlockVariableView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        if constexpr (requires(Index i) { arg.Get()[i]; })
        {
          return arg.Get()[row_in_group];
        }
        else
        {
          return arg.Get();
        }
      }

      /// @brief Get a const element reference for a specific row in this group (Padded Vector-like)
      template<PaddedVectorLike Arg>
      [[gnu::always_inline]]
      decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg[group_ * L + row_in_group];
      }

     public:
      /// @brief Constructor that calculates num_rows_in_group from matrix dimensions
      ConstGroupView(const VectorMatrix& matrix, Index group)
          : matrix_(matrix),
            group_(group)
      {
        // Calculate how many rows are in this group (typically L, except possibly the last group)
        // Optimized: avoid calling NumberOfGroups() which does division+ceil
        Index start_row = group * L;
        num_rows_in_group_ = std::min(L, matrix.x_dim_ - start_row);
      }

      /// @brief Constructor with explicit num_rows_in_group
      ConstGroupView(const VectorMatrix& matrix, Index group, Index num_rows_in_group)
          : matrix_(matrix),
            group_(group),
            num_rows_in_group_(num_rows_in_group)
      {
      }

      /// @brief Returns a grouped const column view whose element base_ pointer is
      ///        precomputed for this ConstGroupView's group.
      GroupedConstColumnView GetConstColumnView(Index column_index) const
      {
        assert(column_index < matrix_.y_dim_ && "column index out of range");
        return { matrix_.data_.data() + (group_ * matrix_.y_dim_ + column_index) * L };
      }

      RowVariable GetRowVariable() const
      {
        // Stack-allocated array of L elements
        return RowVariable();
      }

      /// @brief Assign value to every cell of the caller-owned row-variable temp.
      ///        Handles both flavors of `BlockVariable::Get()`:
      ///          - Array-like (dense RowVariable, and sparse L>1 BlockVariable):
      ///            `Get()` returns `std::array<T, L>&`, so we index it.
      ///          - Scalar (sparse L=1 BlockVariable): `Get()` returns `T&`, so we
      ///            assign directly. Only meaningful when this GroupView's L=1.
      template<BlockVariableView Dst>
      [[gnu::always_inline]]
      void Fill(Dst&& dst, T value) const
      {
        auto& storage = dst.Get();
        if constexpr (requires { storage[Index{ 0 }]; })
        {
          if constexpr (L >= 16)
          {
            std::fill_n(storage.data(), L, value);
          }
          else
          {
            for (Index i = 0; i < L; ++i)
            {
              storage[i] = value;
            }
          }
        }
        else
        {
          static_assert(L == 1, "Scalar BlockVariable::Get() only reachable when L=1");
          storage = value;
        }
      }

      /// @brief Copy src column into the caller-owned row-variable temp.
      ///        See Fill above for the two `Get()` flavors this dispatches over.
      template<BlockVariableView Dst, GroupedDenseMatrixColumnView Src>
      [[gnu::always_inline]]
      void Copy(Dst&& dst, Src&& src) const
      {
        auto& storage = dst.Get();
        if constexpr (requires { storage[Index{ 0 }]; })
        {
          if constexpr (L >= 16)
          {
            std::copy_n(src.base_, L, storage.data());
          }
          else
          {
            for (Index i = 0; i < L; ++i)
            {
              storage[i] = src.base_[i];
            }
          }
        }
        else
        {
          static_assert(L == 1, "Scalar BlockVariable::Get() only reachable when L=1");
          storage = src.base_[0];
        }
      }

      /// @brief Assign value to `vec[group_*L + i]` for every row in this group.
      template<PaddedVectorLike Vec>
      [[gnu::always_inline]]
      void Fill(Vec& vec, T value) const
      {
        const Index start = group_ * L;
        for (Index i = 0; i < L; ++i)
        {
          vec[start + i] = value;
        }
      }

      /// @brief Copy src column into `vec[group_*L .. group_*L + num_rows_in_group_)`.
      template<PaddedVectorLike Vec, GroupedDenseMatrixColumnView Src>
      [[gnu::always_inline]]
      void Copy(Vec& vec, Src&& src) const
      {
        const Index start = group_ * L;
        for (Index i = 0; i < L; ++i)
        {
          vec[start + i] = src.base_[i];
        }
      }

      /// @brief Calls a lambda function for every row in the group (including padded rows)
      template<typename Func, typename... Args>
      void ForEachRow(Func&& func, Args&&... args) const
      {
        for (Index row_in_group = 0; row_in_group < L; ++row_in_group)
        {
          func(GetRowElement(row_in_group, std::forward<Args>(args))...);  // NOLINT(bugprone-use-after-move)
        }
      }

      /// @brief Same as ForEachRow but guaranteed to skip padding rows.
      template<typename Func, typename... Args>
      void ForEachRowStrict(Func&& func, Args&&... args) const
      {
        for (Index row_in_group = 0; row_in_group < num_rows_in_group_; ++row_in_group)
        {
          func(GetRowElement(row_in_group, std::forward<Args>(args))...);  // NOLINT(bugprone-use-after-move)
        }
      }

      /// @brief Apply a reduction to each row in this group. The user's function
      ///        receives its column-view/row-variable arguments plus a trailing
      ///        reference to `reducer.Reference()` as an accumulator, and
      ///        accumulates into it (e.g. `acc += x*x` for a sum, `acc = std::max(acc, x)`
      ///        for a max). Matches ForEachRow's group-iteration shape, including
      ///        operating on padded rows.
      template<typename Reducer, typename Func, typename... Args>
      void Reduce(Reducer reducer, Func&& func, Args&&... args) const
      {
        auto& acc = reducer.Reference();
        for (Index row_in_group = 0; row_in_group < L; ++row_in_group)
        {
          func(GetRowElement(row_in_group, std::forward<Args>(args))..., acc);  // NOLINT(bugprone-use-after-move)
        }
      }

      /// @brief Same as Reduce but guaranteed to skip padding rows.
      template<typename Reducer, typename Func, typename... Args>
      void ReduceStrict(Reducer reducer, Func&& func, Args&&... args) const
      {
        auto& acc = reducer.Reference();
        for (Index row_in_group = 0; row_in_group < num_rows_in_group_; ++row_in_group)
        {
          func(GetRowElement(row_in_group, std::forward<Args>(args))..., acc);  // NOLINT(bugprone-use-after-move)
        }
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

    /// @brief GroupView provides a view of a single group of L rows for iteration
    class GroupView
    {
     public:
      /// @brief Enriched mutable column view returned by GetColumnView on a GroupView.
      ///        See ConstGroupView::GroupedConstColumnView for rationale.
      struct GroupedColumnView
      {
        using category = GroupedDenseMatrixColumnViewTag;
        T* base_;
      };
      /// @brief Const variant, for GetConstColumnView on a mutable GroupView.
      struct GroupedConstColumnView
      {
        using category = GroupedDenseMatrixColumnViewTag;
        const T* base_;
      };

     private:
      VectorMatrix& matrix_;
      Index group_;
      Index num_rows_in_group_;  // May be < L for the last group

      /// @brief Get an element reference for a specific row in this group (ColumnView)
      template<DenseMatrixColumnView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        auto* source_matrix = arg.GetMatrix();
        // VectorMatrix layout: data_[(group * y_dim_ + column) * L + row_in_group]
        return source_matrix->data_[(group_ * source_matrix->y_dim_ + arg.ColumnIndex()) * L + row_in_group];
      }

      /// @brief Get an element reference for a specific row in this group (GroupedColumnView)
      /// Fast path: `base_` already points at row 0 of this group's L-row block.
      template<GroupedDenseMatrixColumnView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg.base_[row_in_group];
      }

      /// @brief Get an element reference for a specific row in this group (BlockVariable).
      ///        See ConstGroupView::GetRowElement above for the subscriptable-vs-scalar
      ///        dispatch rationale.
      template<BlockVariableView Arg>
      [[gnu::always_inline]]
      decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        if constexpr (requires(Index i) { arg.Get()[i]; })
        {
          return arg.Get()[row_in_group];
        }
        else
        {
          return arg.Get();
        }
      }

      /// @brief Get an element reference for a specific row in this group (Vector-like)
      template<PaddedVectorLike Arg>
      [[gnu::always_inline]]
      decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg[group_ * L + row_in_group];
      }

     public:
      /// @brief Constructor that calculates num_rows_in_group from matrix dimensions
      GroupView(VectorMatrix& matrix, Index group)
          : matrix_(matrix),
            group_(group)
      {
        // Calculate how many rows are in this group (typically L, except possibly the last group)
        // Optimized: avoid calling NumberOfGroups() which does division+ceil
        Index start_row = group * L;
        num_rows_in_group_ = std::min(L, matrix.x_dim_ - start_row);
      }

      /// @brief Constructor with explicit num_rows_in_group
      GroupView(VectorMatrix& matrix, Index group, Index num_rows_in_group)
          : matrix_(matrix),
            group_(group),
            num_rows_in_group_(num_rows_in_group)
      {
      }

      operator ConstGroupView() const
      {
        return ConstGroupView(matrix_, group_, num_rows_in_group_);
      }

      /// @brief Returns a grouped const column view whose element base_ pointer is
      ///        precomputed for this GroupView's group.
      GroupedConstColumnView GetConstColumnView(Index column_index) const
      {
        assert(column_index < matrix_.y_dim_ && "column index out of range");
        return { matrix_.data_.data() + (group_ * matrix_.y_dim_ + column_index) * L };
      }

      /// @brief Returns a grouped mutable column view whose element base_ pointer is
      ///        precomputed for this GroupView's group.
      GroupedColumnView GetColumnView(Index column_index) const
      {
        assert(column_index < matrix_.y_dim_ && "column index out of range");
        return { matrix_.data_.data() + (group_ * matrix_.y_dim_ + column_index) * L };
      }

      RowVariable GetRowVariable() const
      {
        // Stack-allocated array of L elements
        return RowVariable();
      }

      [[gnu::always_inline]]
      void Fill(GroupedColumnView view, T value) const
      {
        if constexpr (L >= 16)
        {
          std::fill_n(view.base_, L, value);
        }
        else
        {
          T* dst = view.base_;
          for (Index i = 0; i < L; ++i)
          {
            dst[i] = value;
          }
        }
      }

      template<GroupedDenseMatrixColumnView Src>
      [[gnu::always_inline]]
      void Copy(GroupedColumnView dst_view, Src&& src_view) const
      {
        if constexpr (L >= 16)
        {
          std::copy_n(src_view.base_, L, dst_view.base_);
        }
        else
        {
          T* dst = dst_view.base_;
          const T* src = src_view.base_;
          for (Index i = 0; i < L; ++i)
          {
            dst[i] = src[i];
          }
        }
      }

      /// @brief Copy a per-row vector into dst column within this group.
      template<PaddedVectorLike Src>
      [[gnu::always_inline]]
      void Copy(GroupedColumnView dst_view, Src&& src) const
      {
        T* dst = dst_view.base_;
        const Index start = group_ * L;
        for (Index i = 0; i < num_rows_in_group_; ++i)
        {
          dst[i] = src[start + i];
        }
      }

      /// @brief Assign value to every cell of the caller-owned row-variable temp.
      ///        See ConstGroupView::Fill(Dst&&, T) for the array-vs-scalar
      ///        dispatch rationale.
      template<BlockVariableView Dst>
      [[gnu::always_inline]]
      void Fill(Dst&& dst, T value) const
      {
        auto& storage = dst.Get();
        if constexpr (requires { storage[Index{ 0 }]; })
        {
          if constexpr (L >= 16)
          {
            std::fill_n(storage.data(), L, value);
          }
          else
          {
            for (Index i = 0; i < L; ++i)
            {
              storage[i] = value;
            }
          }
        }
        else
        {
          static_assert(L == 1, "Scalar BlockVariable::Get() only reachable when L=1");
          storage = value;
        }
      }

      /// @brief Copy src column into the caller-owned row-variable temp.
      ///        See ConstGroupView::Copy(Dst&&, Src&&) for the array-vs-scalar
      ///        dispatch rationale.
      template<BlockVariableView Dst, GroupedDenseMatrixColumnView Src>
      [[gnu::always_inline]]
      void Copy(Dst&& dst, Src&& src) const
      {
        auto& storage = dst.Get();
        if constexpr (requires { storage[Index{ 0 }]; })
        {
          if constexpr (L >= 16)
          {
            std::copy_n(src.base_, L, storage.data());
          }
          else
          {
            for (Index i = 0; i < L; ++i)
            {
              storage[i] = src.base_[i];
            }
          }
        }
        else
        {
          static_assert(L == 1, "Scalar BlockVariable::Get() only reachable when L=1");
          storage = src.base_[0];
        }
      }

      /// @brief Assign value to `vec[group_*L + i]` for every real row in this group.
      template<PaddedVectorLike Vec>
      [[gnu::always_inline]]
      void Fill(Vec& vec, T value) const
      {
        const Index start = group_ * L;
        for (Index i = 0; i < L; ++i)
        {
          vec[start + i] = value;
        }
      }

      /// @brief Copy src column into `vec[group_*L .. group_*L + num_rows_in_group_)`.
      template<PaddedVectorLike Vec, GroupedDenseMatrixColumnView Src>
      [[gnu::always_inline]]
      void Copy(Vec& vec, Src&& src) const
      {
        const Index start = group_ * L;
        for (Index i = 0; i < L; ++i)
        {
          vec[start + i] = src.base_[i];
        }
      }

      template<typename Func, typename... Args>
      void ForEachRow(Func&& func, Args&&... args) const
      {
        for (Index row_in_group = 0; row_in_group < L; ++row_in_group)
        {
          func(GetRowElement(row_in_group, std::forward<Args>(args))...);  // NOLINT(bugprone-use-after-move)
        }
      }

      /// @brief Same as ForEachRow but guaranteed to skip padding rows.
      ///        See ConstGroupView::ForEachRowStrict for details.
      template<typename Func, typename... Args>
      void ForEachRowStrict(Func&& func, Args&&... args) const
      {
        for (Index row_in_group = 0; row_in_group < num_rows_in_group_; ++row_in_group)
        {
          func(GetRowElement(row_in_group, std::forward<Args>(args))...);  // NOLINT(bugprone-use-after-move)
        }
      }

      /// @brief Apply a reduction to each row in this group. See
      ///        ConstGroupView::Reduce for details.
      template<typename Reducer, typename Func, typename... Args>
      void Reduce(Reducer reducer, Func&& func, Args&&... args) const
      {
        auto& acc = reducer.Reference();
        for (Index row_in_group = 0; row_in_group < L; ++row_in_group)
        {
          func(GetRowElement(row_in_group, std::forward<Args>(args))..., acc);  // NOLINT(bugprone-use-after-move)
        }
      }

      /// @brief Same as Reduce but guaranteed to skip padding rows.
      template<typename Reducer, typename Func, typename... Args>
      void ReduceStrict(Reducer reducer, Func&& func, Args&&... args) const
      {
        auto& acc = reducer.Reference();
        for (Index row_in_group = 0; row_in_group < num_rows_in_group_; ++row_in_group)
        {
          func(GetRowElement(row_in_group, std::forward<Args>(args))..., acc);  // NOLINT(bugprone-use-after-move)
        }
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

    /// @brief Create a function that can be applied to vector matrices and vectors
    ///
    /// Creates a reusable callable that validates matrix dimensions and applies a user function
    /// across row groups. The function iterates over groups of L rows at a time for vectorization,
    /// where L is the compile-time template parameter.
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
    /// @throws std::system_error if column counts don't match at creation, or if at invocation time:
    ///         matrices/vectors have mismatched row counts, column counts don't match creation,
    ///         or dimensions mismatch
    ///
    /// @tparam UseView When true (default), vector args are converted to their View/ConstView
    ///                 via `arg.GetView()` before being handed to the lambda so lambdas whose
    ///                 parameter is declared `Vector::ViewType`/`Vector::ConstViewType` see a
    ///                 lightweight view (mirrors the Kokkos MakeHandle path).  When false,
    ///                 vector args are passed through unchanged; use this for HostFunction
    ///                 where the arg may be a KokkosPaddedVector whose GetView() returns a
    ///                 device view unusable on host.  Both host PaddedVector and
    ///                 KokkosPaddedVector satisfy PaddedVectorLike (operator[], size,
    ///                 PaddedSize), so GroupView::GetRowElement handles both.
    template<bool UseView = true, typename Func, typename... Args>
    static auto Function(Func&& func, Args&... args)
    {
      // Capture column counts for matrices at creation time using helper
      // Row counts can differ between args at creation, but must match at invocation
      auto populate_cols = [](auto&... args_inner)
      {
        std::vector<Index> cols(sizeof...(args_inner));
        Index idx = 0;
        (
            [&](auto& arg)
            {
              using ArgType = std::remove_cvref_t<decltype(arg)>;
              if constexpr (PaddedVectorLike<ArgType>)
              {
                cols[idx] = 0;  // Not used for vectors
              }
              else
              {
                cols[idx] = arg.NumColumns();
              }
              ++idx;
            }(args_inner),
            ...);
        return cols;
      };

      std::vector<Index> num_cols = populate_cols(args...);

      // Store in variable to ensure fold expression completes before lambda construction
      auto result = [func = std::forward<Func>(func), num_cols = std::move(num_cols)](auto&&... invoked_args) mutable
      {
        // Validate dimensions and determine row count in a single pass
        Index num_rows = 0;
        bool found_first = false;
        Index idx = 0;

        (
            [&](auto& arg)
            {
              using ArgType = std::remove_cvref_t<decltype(arg)>;

              if constexpr (PaddedVectorLike<ArgType>)
              {
                // Vector - validate size matches row count
                if (!found_first)
                {
                  num_rows = arg.size();
                  found_first = true;
                }
                else if (arg.size() != num_rows)
                {
                  throw MicmException(
                      MICM_ERROR_CATEGORY_MATRIX,
                      MICM_MATRIX_ERROR_CODE_INVALID_VECTOR,
                      "Vector size must match matrix row count. Expected " + std::to_string(num_rows) +
                          " elements but got " + std::to_string(arg.size()));
                }
              }
              else if constexpr (requires {
                                   arg.NumRows();
                                   arg.NumColumns();
                                 })
              {
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
                    throw MicmException(
                        MICM_ERROR_CATEGORY_MATRIX,
                        MICM_MATRIX_ERROR_CODE_INVALID_VECTOR,
                        "All matrices must have the same number of rows when invoking function. Expected " +
                            std::to_string(num_rows) + " rows but got " + std::to_string(arg.NumRows()));
                  }
                }

                // Always validate column count against captured value
                if (arg.NumColumns() != num_cols[idx])
                {
                  throw MicmException(
                      MICM_ERROR_CATEGORY_MATRIX,
                      MICM_MATRIX_ERROR_CODE_INVALID_VECTOR,
                      "Matrix column count does not match. Expected " + std::to_string(num_cols[idx]) + " columns but got " +
                          std::to_string(arg.NumColumns()));
                }
              }
              ++idx;
            }(invoked_args),
            ...);

        // Iterate over groups, processing L rows at a time
        Index num_complete_groups = (num_rows / L);
        for (Index group = 0; group < num_complete_groups; ++group)
        {
          // Use ConstGroupView if matrix is const, otherwise use GroupView
          // For vectors, either take a View/ConstView (UseView=true, matches lambda param
          // typed as Vector::ViewType/ConstViewType) or pass through raw (UseView=false,
          // required for HostFunction which may receive KokkosPaddedVector arguments whose
          // GetView() returns a device view).
          func(
              [&](auto&& arg) -> decltype(auto)
              {
                using ArgType = std::remove_reference_t<decltype(arg)>;
                using ArgTypeNoConst = std::remove_const_t<ArgType>;
                if constexpr (PaddedVectorLike<std::remove_cvref_t<ArgType>>)
                {
                  if constexpr (UseView)
                  {
                    return std::forward<decltype(arg)>(arg).GetView();
                  }
                  else
                  {
                    return (std::forward<decltype(arg)>(arg));
                  }
                }
                else
                {
                  // Matrix: create appropriate GroupView
                  if constexpr (std::is_const_v<ArgType>)
                  {
                    return typename ArgTypeNoConst::ConstHostGroupView(arg, group, L);
                  }
                  else
                  {
                    return typename ArgTypeNoConst::HostGroupView(arg, group, L);
                  }
                }
              }(invoked_args)...);
        }

        // Process remaining rows (if num_rows is not a multiple of L)
        Index remaining = num_rows % L;
        if (remaining > 0)
        {
          // Use ConstGroupView if matrix is const, otherwise use GroupView
          // For vectors, see the matching complete-group case for UseView vs pass-through
          // rationale.
          func(
              [&](auto&& arg) -> decltype(auto)
              {
                using ArgType = std::remove_reference_t<decltype(arg)>;
                using ArgTypeNoConst = std::remove_const_t<ArgType>;
                if constexpr (PaddedVectorLike<std::remove_cvref_t<ArgType>>)
                {
                  if constexpr (UseView)
                  {
                    return std::forward<decltype(arg)>(arg).GetView();
                  }
                  else
                  {
                    return (std::forward<decltype(arg)>(arg));
                  }
                }
                else
                {
                  // Matrix: create appropriate GroupView
                  if constexpr (std::is_const_v<ArgType>)
                  {
                    return typename ArgTypeNoConst::ConstHostGroupView(arg, num_complete_groups, remaining);
                  }
                  else
                  {
                    return typename ArgTypeNoConst::HostGroupView(arg, num_complete_groups, remaining);
                  }
                }
              }(invoked_args)...);
        }
      };
      return result;
    }

    template<typename Func, typename... Args>
    static auto HostFunction(Func&& func, Args&... args)
    {
      // UseView=false: pass vector args through unchanged.  A KokkosPaddedVector's
      // GetView() returns a device view unusable on host; passing the padded vector
      // itself lets GroupView::GetRowElement dispatch to the PaddedVectorLike overload
      // which just indexes via operator[].
      return Function<false>(std::forward<Func>(func), args...);
    }

   private:
    /// @brief Get an element reference for a row (ColumnView)
    template<DenseMatrixColumnView Arg>
    [[gnu::always_inline]]
    decltype(auto) GetRowElement(Index row, Index group, Index row_in_group, Arg&& arg)
    {
      auto* source_matrix = arg.GetMatrix();
      // VectorMatrix layout: data_[(group * y_dim_ + column) * L + row_in_group]
      return source_matrix->data_[(group * source_matrix->y_dim_ + arg.ColumnIndex()) * L + row_in_group];
    }

    /// @brief Get an element reference for a row (BlockVariable).
    ///        See GroupView::GetRowElement above for the subscriptable-vs-scalar
    ///        dispatch rationale.
    template<BlockVariableView Arg>
    [[gnu::always_inline]]
    decltype(auto) GetRowElement(Index /*row*/, Index /*group*/, Index row_in_group, Arg&& arg)
    {
      if constexpr (requires(Index i) { arg.Get()[i]; })
      {
        return arg.Get()[row_in_group];
      }
      else
      {
        return arg.Get();
      }
    }

    /// @brief Get an element reference for a row (Vector-like)
    template<PaddedVectorLike Arg>
    [[gnu::always_inline]]
    decltype(auto) GetRowElement(Index row, Index group, Index row_in_group, Arg&& arg)
    {
      return arg[row];
    }

    /// @brief Get a const element reference for a row (ColumnView) - const version
    template<DenseMatrixColumnView Arg>
    [[gnu::always_inline]]
    decltype(auto) GetRowElement(Index row, Index group, Index row_in_group, Arg&& arg) const
    {
      auto* source_matrix = arg.GetMatrix();
      // VectorMatrix layout: data_[(group * y_dim_ + column) * L + row_in_group]
      return source_matrix->data_[(group * source_matrix->y_dim_ + arg.ColumnIndex()) * L + row_in_group];
    }

    /// @brief Get a const element reference for a row (BlockVariable) - const version.
    ///        See GroupView::GetRowElement above for the subscriptable-vs-scalar
    ///        dispatch rationale.
    template<BlockVariableView Arg>
    [[gnu::always_inline]]
    decltype(auto) GetRowElement(Index /*row*/, Index /*group*/, Index row_in_group, Arg&& arg) const
    {
      if constexpr (requires(Index i) { arg.Get()[i]; })
      {
        return arg.Get()[row_in_group];
      }
      else
      {
        return arg.Get();
      }
    }

    /// @brief Get a const element reference for a row (Vector-like) - const version
    template<PaddedVectorLike Arg>
    [[gnu::always_inline]]
    decltype(auto) GetRowElement(Index row, Index group, Index row_in_group, Arg&& arg) const
    {
      return arg[row];
    }
  };

  // ============================================================================
  // Grouping Strategy Specialization
  // ============================================================================

  /// @brief VectorMatrix uses simple grouping when L==1, tiered grouping when L>1
  template<typename T, Index L>
  struct GroupingStrategy<VectorMatrix<T, L>>
  {
    using type = std::conditional_t<L == 1, SimpleGroupingTag, TieredGroupingTag>;
  };

}  // namespace micm
