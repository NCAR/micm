// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/kokkos/util/kokkos_views.hpp>
#include <micm/util/vector_matrix.hpp>

#include <Kokkos_Core.hpp>
#include <vector>

namespace micm
{
  /// @brief Provides a Kokkos implementation to the VectorMatrix functionality.
  ///
  /// Inherits from VectorMatrix (the MICM host-side data layout) and maintains
  /// a Kokkos::View as a device-side mirror. The caller must explicitly call
  /// CopyToDevice() / CopyToHost() to synchronize, matching the CUDA matrix pattern.
  template<class T, Index L = MICM_DEFAULT_VECTOR_SIZE>
  class KokkosDenseMatrix : public VectorMatrix<T, L>
  {
   public:
    static constexpr Index GroupVectorSize()
    {
      return L;
    }
    using value_type = T;
    using ViewType = Kokkos::View<T*>;
    using HostViewType = typename ViewType::host_mirror_type;
    using TeamPolicyType = Kokkos::TeamPolicy<>;
    using TeamMember = typename TeamPolicyType::member_type;

   private:
    /// Device-side (or unified) view — the Kokkos mirror of VectorMatrix::data_
    ViewType view_;

    /// @brief Device-safe handle for a mutable KokkosDenseMatrix argument to
    ///        Function()/ForEachRow().
    ///
    /// Only the flat Kokkos::View and column count are captured -- both are trivially
    /// copyable -- so this handle (rather than the matrix object itself, which owns a
    /// non-trivially-copyable host std::vector) is what gets captured by value into a
    /// device lambda. Private to KokkosDenseMatrix -- only ever used internally by
    /// MakeHandle()/BuildGroupView() below.
    struct DenseMatrixHandle
    {
      Kokkos::View<T*> view;
      Index y_dim;
    };

    /// @brief Const variant of DenseMatrixHandle. See DenseMatrixHandle for details.
    struct ConstDenseMatrixHandle
    {
      Kokkos::View<const T*> view;
      Index y_dim;
    };

    /// @brief Build a device-safe handle for one Function()/ForEachRow() argument.
    ///
    /// Kokkos matrix arguments are reduced to their (trivially copyable) View + column
    /// count. VectorLike arguments (e.g. a host std::vector) are forwarded unchanged --
    /// this only works correctly on host-accessible backends (Serial, OpenMP), and even
    /// there, only for read-only usage: KOKKOS_LAMBDA captures everything by value
    /// (`[=]`), and per the C++ standard, capturing a *reference*-typed local by value
    /// copies the *referenced object*, so writes into a VectorLike argument inside the
    /// kernel are NOT observable by the caller afterward. This is a known limitation,
    /// deferred to a follow-up that replaces VectorLike arguments with a typed
    /// device-aware vector view (see the plan's "Solver lambdas on GPU / VectorLike
    /// args" design note).
    template<typename Arg>
    static auto MakeHandle(Arg&& arg)
    {
      using ArgType = std::remove_reference_t<Arg>;
      if constexpr (VectorLike<std::remove_cvref_t<ArgType>>)
      {
        return std::forward<Arg>(arg);
      }
      else if constexpr (std::is_const_v<ArgType>)
      {
        return ConstDenseMatrixHandle{ arg.GetView(), arg.NumColumns() };
      }
      else
      {
        return DenseMatrixHandle{ arg.GetView(), arg.NumColumns() };
      }
    }

    /// @brief Construct the appropriate GroupView/ConstGroupView (or forward a
    ///        VectorLike argument unchanged) for one handle produced by
    ///        MakeHandle(). Runs on-device (called from within a KOKKOS_LAMBDA).
    template<typename Handle>
    KOKKOS_INLINE_FUNCTION static decltype(auto) BuildGroupView(
        Handle&& handle,
        Index group,
        Index count,
        const TeamMember& team)
    {
      using HandleType = std::remove_cvref_t<Handle>;
      if constexpr (std::is_same_v<HandleType, DenseMatrixHandle>)
      {
        return GroupView(handle.view, group, handle.y_dim, count, team);
      }
      else if constexpr (std::is_same_v<HandleType, ConstDenseMatrixHandle>)
      {
        return ConstGroupView(handle.view, group, handle.y_dim, count, team);
      }
      else
      {
        // VectorLike: forward through as-is (see MakeHandle() note)
        return std::forward<Handle>(handle);
      }
    }

   public:
    KokkosDenseMatrix()
        : VectorMatrix<T, L>()
    {
    }

    KokkosDenseMatrix(Index x_dim, Index y_dim)
        : VectorMatrix<T, L>(x_dim, y_dim),
          view_("dense_matrix", VectorMatrix<T, L>(x_dim, y_dim).AsVector().size())
    {
    }

    KokkosDenseMatrix(Index x_dim, Index y_dim, T initial_value)
        : VectorMatrix<T, L>(x_dim, y_dim, initial_value),
          view_("dense_matrix", VectorMatrix<T, L>(x_dim, y_dim).AsVector().size())
    {
      Kokkos::deep_copy(view_, initial_value);
    }

    KokkosDenseMatrix(const std::vector<std::vector<T>>& other)
        : VectorMatrix<T, L>(other),
          view_("dense_matrix", this->data_.size())
    {
      CopyToDevice();
    }

    /// @brief Copy host data (MICM's data_) to the device view
    void CopyToDevice()
    {
      if (view_.extent(0) != this->data_.size())
      {
        view_ = ViewType("dense_matrix", this->data_.size());
      }
      auto h_view = Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
          this->data_.data(), this->data_.size());
      Kokkos::deep_copy(view_, h_view);
    }

    /// @brief Copy device view data back to host (MICM's data_)
    void CopyToHost()
    {
      if (view_.extent(0) != 0)
      {
        auto h_view = Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            this->data_.data(), this->data_.size());
        Kokkos::deep_copy(h_view, view_);
      }
    }

    ViewType GetView() const
    {
      return view_;
    }

    /// @brief Set every element on the device to a given value
    void Fill(T val)
    {
      if (view_.extent(0) != this->data_.size())
      {
        view_ = ViewType("dense_matrix", this->data_.size());
      }
      Kokkos::deep_copy(view_, val);
    }

    /// @brief For each element in the KokkosDenseMatrix x and y, perform y = alpha * x + y,
    ///        where alpha is a scalar constant. Runs on-device.
    ///
    /// Only touches the matrix's real (non-padding) cells.
    /// @param alpha The scaling scalar to apply to the KokkosDenseMatrix x
    /// @param x The input KokkosDenseMatrix
    void Axpy(const Real& alpha, const KokkosDenseMatrix& x)
    {
      ViewType y_view = view_;
      ViewType x_view = x.view_;
      const Index y_dim = this->NumColumns();
      const Index n = static_cast<Index>(std::floor(this->NumRows() / (double)L)) * L * y_dim;
      Kokkos::parallel_for(
          "KokkosDenseMatrix::Axpy",
          Kokkos::RangePolicy<>(0, n),
          KOKKOS_LAMBDA(const Index i) { y_view(i) += alpha * x_view(i); });
      const Index l = this->NumRows() % L;
      if (l > 0)
      {
        Kokkos::parallel_for(
            "KokkosDenseMatrix::Axpy(tail)",
            Kokkos::RangePolicy<>(0, y_dim * l),
            KOKKOS_LAMBDA(const Index idx)
            {
              const Index flat = n + (idx / l) * L + (idx % l);
              y_view(flat) += alpha * x_view(flat);
            });
      }
    }

    /// @brief For each element of the matrix, perform y = max(y, x), where x is a scalar constant.
    ///
    /// Touches every stored cell, including any trailing padding cells.
    void Max(const T& x)
    {
      ViewType y_view = view_;
      Kokkos::parallel_for(
          "KokkosDenseMatrix::Max",
          Kokkos::RangePolicy<>(0, y_view.extent(0)),
          KOKKOS_LAMBDA(const Index i) { y_view(i) = y_view(i) > x ? y_view(i) : x; });
    }

    /// @brief For each element of the matrix, perform y = min(y, x), where x is a scalar constant.
    ///
    /// Touches every stored cell, including any trailing padding cells.
    void Min(const T& x)
    {
      ViewType y_view = view_;
      Kokkos::parallel_for(
          "KokkosDenseMatrix::Min",
          Kokkos::RangePolicy<>(0, y_view.extent(0)),
          KOKKOS_LAMBDA(const Index i) { y_view(i) = y_view(i) < x ? y_view(i) : x; });
    }

    /// @brief Copy the device data from the other Kokkos dense matrix into this one
    void Copy(const KokkosDenseMatrix& other)
    {
      if (other.view_.extent(0) != view_.extent(0))
      {
        throw std::runtime_error("Both Kokkos dense matrices must have the same size.");
      }
      Kokkos::deep_copy(view_, other.view_);
    }

    /// @brief Swap the device data from the other Kokkos dense matrix into this one.
    void Swap(KokkosDenseMatrix& other)
    {
      if (other.view_.extent(0) != view_.extent(0))
      {
        throw std::runtime_error("Both Kokkos dense matrices must have the same size.");
      }
      std::swap(view_, other.view_);
      this->data_.swap(other.data_);
    }

    KokkosDenseMatrix& operator=(T val)
    {
      Fill(val);
      return *this;
    }

    /// @brief Apply a two-argument element-wise function on-device.
    ///
    /// Only touches real (non-padding) cells.
    template<typename Func>
    void ForEach(Func&& f, const KokkosDenseMatrix& a)
    {
      ViewType y_view = view_;
      ViewType a_view = a.view_;
      const Index y_dim = this->NumColumns();
      const Index n = static_cast<Index>(std::floor(this->NumRows() / (double)L)) * L * y_dim;
      Kokkos::parallel_for(
          "KokkosDenseMatrix::ForEach",
          Kokkos::RangePolicy<>(0, n),
          KOKKOS_LAMBDA(const Index i) { f(y_view(i), a_view(i)); });
      const Index l = this->NumRows() % L;
      if (l > 0)
      {
        Kokkos::parallel_for(
            "KokkosDenseMatrix::ForEach(tail)",
            Kokkos::RangePolicy<>(0, y_dim * l),
            KOKKOS_LAMBDA(const Index idx)
            {
              const Index flat = n + (idx / l) * L + (idx % l);
              f(y_view(flat), a_view(flat));
            });
      }
    }

    /// @brief Apply a three-argument element-wise function on-device.
    ///
    /// Only touches real (non-padding) cells; see Axpy() note.
    template<typename Func>
    void ForEach(Func&& f, const KokkosDenseMatrix& a, const KokkosDenseMatrix& b)
    {
      ViewType y_view = view_;
      ViewType a_view = a.view_;
      ViewType b_view = b.view_;
      const Index y_dim = this->NumColumns();
      const Index n = static_cast<Index>(std::floor(this->NumRows() / (double)L)) * L * y_dim;
      Kokkos::parallel_for(
          "KokkosDenseMatrix::ForEach",
          Kokkos::RangePolicy<>(0, n),
          KOKKOS_LAMBDA(const Index i) { f(y_view(i), a_view(i), b_view(i)); });
      const Index l = this->NumRows() % L;
      if (l > 0)
      {
        Kokkos::parallel_for(
            "KokkosDenseMatrix::ForEach(tail)",
            Kokkos::RangePolicy<>(0, y_dim * l),
            KOKKOS_LAMBDA(const Index idx)
            {
              const Index flat = n + (idx / l) * L + (idx % l);
              f(y_view(flat), a_view(flat), b_view(flat));
            });
      }
    }

    /// @brief ConstGroupView provides a team-parallel const view of a single group of L
    ///        rows for iteration on-device.
    class ConstGroupView
    {
     public:
      using GroupedConstColumnView = micm::KokkosGroupedConstColumnView<T>;

     private:
      Kokkos::View<const T*> view_;
      Index group_;
      Index y_dim_;
      Index num_rows_in_group_;
      TeamMember team_;

      template<DenseMatrixColumnView Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg.Data()[(group_ * arg.YDim() + arg.ColumnIndex()) * L + row_in_group];
      }

      template<GroupedDenseMatrixColumnView Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg.base_[row_in_group];
      }

      template<BlockVariableView Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg.Get()[row_in_group];
      }

      template<VectorLike Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg[group_ * L + row_in_group];
      }

     public:
      KOKKOS_INLINE_FUNCTION
      ConstGroupView(Kokkos::View<const T*> view, Index group, Index y_dim, Index num_rows_in_group, const TeamMember& team)
          : view_(view),
            group_(group),
            y_dim_(y_dim),
            num_rows_in_group_(num_rows_in_group),
            team_(team)
      {
      }

      KOKKOS_INLINE_FUNCTION GroupedConstColumnView GetConstColumnView(Index column_index) const
      {
        return { view_.data() + (group_ * y_dim_ + column_index) * L };
      }

      KOKKOS_INLINE_FUNCTION KokkosRowVariable<T, L> GetRowVariable() const
      {
        return KokkosRowVariable<T, L>();
      }

      /// @brief Assign value to dst.
      template<BlockVariableView Dst>
      KOKKOS_INLINE_FUNCTION void Fill(Dst&& dst, T value) const
      {
        auto& storage = dst.Get();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_, L), [&](const Index i) { storage[i] = value; });
        team_.team_barrier();
      }

      /// @brief Copy src into dst.
      template<BlockVariableView Dst, GroupedDenseMatrixColumnView Src>
      KOKKOS_INLINE_FUNCTION void Copy(Dst&& dst, Src&& src) const
      {
        auto& storage = dst.Get();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_, L), [&](const Index i) { storage[i] = src.base_[i]; });
        team_.team_barrier();
      }

      /// @brief Assign value to every element in the vector.
      template<VectorLike Vec>
      KOKKOS_INLINE_FUNCTION void Fill(Vec& vec, T value) const
      {
        const Index start = group_ * L;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_rows_in_group_), [&](const Index i) { vec[start + i] = value; });
        team_.team_barrier();
      }

      /// @brief Copy src into vec.
      template<VectorLike Vec, GroupedDenseMatrixColumnView Src>
      KOKKOS_INLINE_FUNCTION void Copy(Vec& vec, Src&& src) const
      {
        const Index start = group_ * L;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_rows_in_group_), [&](const Index i) { vec[start + i] = src.base_[i]; });
        team_.team_barrier();
      }

      /// @brief Apply the provided function to every row in the matrix.
      template<typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void ForEachRow(Func&& func, Args&&... args) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, L),
            [&](const Index row_in_group) { func(GetRowElement(row_in_group, args)...); });
        team_.team_barrier();
      }

      /// @brief Same as ForEachRow but guaranteed to skip padding rows.
      template<typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void ForEachRowStrict(Func&& func, Args&&... args) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_rows_in_group_),
            [&](const Index row_in_group) { func(GetRowElement(row_in_group, args)...); });
        team_.team_barrier();
      }
    };

    /// @brief GroupView provides a team-parallel view of a single group of L rows for
    ///        iteration on-device.
    class GroupView
    {
     public:
      using GroupedColumnView = micm::KokkosGroupedColumnView<T>;
      using GroupedConstColumnView = micm::KokkosGroupedConstColumnView<T>;

     private:
      ViewType view_;
      Index group_;
      Index y_dim_;
      Index num_rows_in_group_;
      TeamMember team_;

      template<DenseMatrixColumnView Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg.Data()[(group_ * arg.YDim() + arg.ColumnIndex()) * L + row_in_group];
      }

      template<GroupedDenseMatrixColumnView Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg.base_[row_in_group];
      }

      template<BlockVariableView Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg.Get()[row_in_group];
      }

      template<VectorLike Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg[group_ * L + row_in_group];
      }

     public:
      KOKKOS_INLINE_FUNCTION
      GroupView(ViewType view, Index group, Index y_dim, Index num_rows_in_group, const TeamMember& team)
          : view_(view),
            group_(group),
            y_dim_(y_dim),
            num_rows_in_group_(num_rows_in_group),
            team_(team)
      {
      }

      KOKKOS_INLINE_FUNCTION GroupedConstColumnView GetConstColumnView(Index column_index) const
      {
        return { view_.data() + (group_ * y_dim_ + column_index) * L };
      }

      KOKKOS_INLINE_FUNCTION GroupedColumnView GetColumnView(Index column_index) const
      {
        return { view_.data() + (group_ * y_dim_ + column_index) * L };
      }

      KOKKOS_INLINE_FUNCTION KokkosRowVariable<T, L> GetRowVariable() const
      {
        return KokkosRowVariable<T, L>();
      }

      /// @brief Assign value to view.
      KOKKOS_INLINE_FUNCTION void Fill(GroupedColumnView view, T value) const
      {
        T* dst = view.base_;
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_, L), [&](const Index i) { dst[i] = value; });
        team_.team_barrier();
      }

      /// @brief Copy src column into dst_view.
      template<GroupedDenseMatrixColumnView Src>
      KOKKOS_INLINE_FUNCTION void Copy(GroupedColumnView dst_view, Src&& src_view) const
      {
        T* dst = dst_view.base_;
        const T* src = src_view.base_;
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_, L), [&](const Index i) { dst[i] = src[i]; });
        team_.team_barrier();
      }

      /// @brief Copy src into dst_view.
      template<VectorLike Src>
      KOKKOS_INLINE_FUNCTION void Copy(GroupedColumnView dst_view, Src&& src) const
      {
        T* dst = dst_view.base_;
        const Index start = group_ * L;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_rows_in_group_), [&](const Index i) { dst[i] = src[start + i]; });
        team_.team_barrier();
      }

      /// @brief Assign value to dst.
      template<BlockVariableView Dst>
      KOKKOS_INLINE_FUNCTION void Fill(Dst&& dst, T value) const
      {
        auto& storage = dst.Get();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_, L), [&](const Index i) { storage[i] = value; });
        team_.team_barrier();
      }

      /// @brief Copy src into dst.
      template<BlockVariableView Dst, GroupedDenseMatrixColumnView Src>
      KOKKOS_INLINE_FUNCTION void Copy(Dst&& dst, Src&& src) const
      {
        auto& storage = dst.Get();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_, L), [&](const Index i) { storage[i] = src.base_[i]; });
        team_.team_barrier();
      }

      /// @brief Assign value to all vec elements.
      template<VectorLike Vec>
      KOKKOS_INLINE_FUNCTION void Fill(Vec& vec, T value) const
      {
        const Index start = group_ * L;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_rows_in_group_), [&](const Index i) { vec[start + i] = value; });
        team_.team_barrier();
      }

      /// @brief Copy src into vec.
      template<VectorLike Vec, GroupedDenseMatrixColumnView Src>
      KOKKOS_INLINE_FUNCTION void Copy(Vec& vec, Src&& src) const
      {
        const Index start = group_ * L;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_rows_in_group_), [&](const Index i) { vec[start + i] = src.base_[i]; });
        team_.team_barrier();
      }

      /// @brief Apply the provided function to every row in the matrix.
      template<typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void ForEachRow(Func&& func, Args&&... args) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, L),
            [&](const Index row_in_group) { func(GetRowElement(row_in_group, args)...); });
        team_.team_barrier();
      }

      /// @brief Same as ForEachRow but guaranteed to skip padding rows.
      template<typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void ForEachRowStrict(Func&& func, Args&&... args) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_rows_in_group_),
            [&](const Index row_in_group) { func(GetRowElement(row_in_group, args)...); });
        team_.team_barrier();
      }
    };

    /// @brief Create a function that can be applied to Kokkos dense matrices and
    ///        vectors, executing on-device using team parallelism.
    template<typename Func, typename... Args>
    static auto Function(Func&& func, Args&... args)
    {
      auto populate_cols = [](auto&... args_inner)
      {
        std::vector<Index> cols(sizeof...(args_inner));
        Index idx = 0;
        (
            [&](auto& arg)
            {
              using ArgType = std::remove_cvref_t<decltype(arg)>;
              if constexpr (VectorLike<ArgType>)
              {
                cols[idx] = 0;
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

      auto result = [func = std::forward<Func>(func), num_cols = std::move(num_cols)](auto&&... invoked_args) mutable
      {
        Index num_rows = 0;
        bool found_first = false;
        Index idx = 0;

        (
            [&](auto& arg)
            {
              using ArgType = std::remove_cvref_t<decltype(arg)>;

              if constexpr (VectorLike<ArgType>)
              {
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
                if (!found_first)
                {
                  num_rows = arg.NumRows();
                  found_first = true;
                }
                else if (arg.NumRows() != num_rows)
                {
                  throw MicmException(
                      MICM_ERROR_CATEGORY_MATRIX,
                      MICM_MATRIX_ERROR_CODE_INVALID_VECTOR,
                      "All matrices must have the same number of rows when invoking function. Expected " +
                          std::to_string(num_rows) + " rows but got " + std::to_string(arg.NumRows()));
                }

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

        Index num_complete_groups = static_cast<Index>(std::floor(num_rows / (double)L));
        Index remaining = num_rows % L;

        // Reduce each argument to a device-safe handle (Kokkos::View + column count for
        // matrices; forwarded as-is for VectorLike) *before* entering device code, then
        // hand the handles to a templated lambda so the pack stays a set of plain,
        // individually-typed parameters rather than a (device-unfriendly) std::tuple.
        //
        // Always dispatch via Kokkos::TeamPolicy (never a bare RangePolicy), even when
        // L == 1: GroupView/ConstGroupView require a real TeamMember (their
        // ForEachRow/Fill/Copy helpers use Kokkos::TeamThreadRange(team_, ...) and
        // team_.team_barrier()), and Kokkos's TeamMember types are not
        // default-constructible, so there is no way to hand GroupView a placeholder
        // "no team" value. For L == 1, num_complete_groups == num_rows and
        // remaining == 0, so this reduces to one team of size 1 per row -- correct,
        // if not maximally GPU-efficient (see the plan's team-size note).
        [&]<typename... Handles>(Handles&&... handles)
        {
          if (num_complete_groups > 0)
          {
            TeamPolicyType policy(static_cast<int>(num_complete_groups), Kokkos::AUTO);
            Kokkos::parallel_for(
                "KokkosDenseMatrix::Function",
                policy,
                KOKKOS_LAMBDA(const TeamMember& team)
                {
                  const Index group = static_cast<Index>(team.league_rank());
                  func(BuildGroupView(handles, group, L, team)...);
                });
          }
          if (remaining > 0)
          {
            TeamPolicyType tail_policy(1, Kokkos::AUTO);
            Kokkos::parallel_for(
                "KokkosDenseMatrix::Function(tail)",
                tail_policy,
                KOKKOS_LAMBDA(const TeamMember& team)
                { func(BuildGroupView(handles, num_complete_groups, remaining, team)...); });
          }
        }(MakeHandle(invoked_args)...);
      };
      return result;
    }

    /// @brief Apply a function to each row of the matrix, executing on-device using
    ///        team parallelism (one team per row-group of L rows).
    /// @tparam Func The lambda/function type
    /// @tparam Args The types of the column view / vector arguments
    /// @param func The function to apply to each row
    /// @param args Column views, row variables, or vectors
    template<typename Func, typename... Args>
    void ForEachRow(Func&& func, Args&&... args)
    {
      const Index num_rows = this->NumRows();
      const Index y_dim = this->NumColumns();
      ViewType view = view_;

      if constexpr (L == 1)
      {
        if (num_rows > 0)
        {
          Kokkos::parallel_for(
              "KokkosDenseMatrix::ForEachRow",
              Kokkos::RangePolicy<>(0, num_rows),
              KOKKOS_LAMBDA(const Index row) { func(GetTopLevelRowElement(view, y_dim, row, args)...); });
        }
        return;
      }

      const Index num_complete_groups = static_cast<Index>(std::floor(num_rows / (double)L));
      const Index remaining = num_rows % L;

      if (num_complete_groups > 0)
      {
        TeamPolicyType policy(static_cast<int>(num_complete_groups), Kokkos::AUTO);
        Kokkos::parallel_for(
            "KokkosDenseMatrix::ForEachRow",
            policy,
            KOKKOS_LAMBDA(const TeamMember& team)
            {
              const Index group = static_cast<Index>(team.league_rank());
              Kokkos::parallel_for(
                  Kokkos::TeamThreadRange(team, L),
                  [&](const Index row_in_group)
                  {
                    const Index row = group * L + row_in_group;
                    func(GetTopLevelRowElement(view, y_dim, row, args)...);
                  });
            });
      }
      if (remaining > 0)
      {
        TeamPolicyType tail_policy(1, Kokkos::AUTO);
        Kokkos::parallel_for(
            "KokkosDenseMatrix::ForEachRow(tail)",
            tail_policy,
            KOKKOS_LAMBDA(const TeamMember& team)
            {
              Kokkos::parallel_for(
                  Kokkos::TeamThreadRange(team, remaining),
                  [&](const Index row_in_group)
                  {
                    const Index row = num_complete_groups * L + row_in_group;
                    func(GetTopLevelRowElement(view, y_dim, row, args)...);
                  });
            });
      }
    }

   private:
    /// @brief Get an element reference for a row at the (ungrouped) matrix level.
    ///        Used by the matrix-level ForEachRow() override.
    template<DenseMatrixColumnView Arg>
    KOKKOS_INLINE_FUNCTION static decltype(auto) GetTopLevelRowElement(ViewType, Index, Index row, Arg&& arg)
    {
      return arg.Data()[(row / L * arg.YDim() + arg.ColumnIndex()) * L + row % L];
    }

    template<BlockVariableView Arg>
    KOKKOS_INLINE_FUNCTION static decltype(auto) GetTopLevelRowElement(ViewType, Index, Index row, Arg&& arg)
    {
      return arg.Get()[row % L];
    }

    template<VectorLike Arg>
    KOKKOS_INLINE_FUNCTION static decltype(auto) GetTopLevelRowElement(ViewType, Index, Index row, Arg&& arg)
    {
      return arg[row];
    }
  };
}  // namespace micm
