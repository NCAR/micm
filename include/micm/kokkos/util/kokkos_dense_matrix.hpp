// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/kokkos/util/kokkos_padded_vector.hpp>
#include <micm/kokkos/util/kokkos_view_category.hpp>
#include <micm/kokkos/util/kokkos_views.hpp>
#include <micm/util/types.hpp>
#include <micm/util/reducers.hpp>
#include <micm/util/vector_matrix.hpp>

#include <Kokkos_Core.hpp>
#include <vector>

namespace micm
{
  namespace detail
  {
    /// @brief Maps a micm reducer type (Sum/Max/LOr/LAnd) to the matching Kokkos
    ///        reducer type used to actually drive Kokkos::parallel_reduce.
    template<typename Reducer>
    struct ToKokkosReducer;

    template<typename T>
    struct ToKokkosReducer<micm::Sum<T>>
    {
      using type = Kokkos::Sum<T>;
    };

    template<typename T>
    struct ToKokkosReducer<micm::Max<T>>
    {
      using type = Kokkos::Max<T>;
    };

    template<>
    struct ToKokkosReducer<micm::LOr>
    {
      using type = Kokkos::LOr<bool>;
    };

    template<>
    struct ToKokkosReducer<micm::LAnd>
    {
      using type = Kokkos::LAnd<bool>;
    };

    // ----------------------------------------------------------------
    // Device-compatible tuple for KOKKOS_LAMBDA-safe variadic dispatch
    // ----------------------------------------------------------------
    // EBO-based heterogeneous tuple: no STL runtime dependencies,
    // works inside KOKKOS_FUNCTION contexts.
    // Used to bundle Function()/ForEachRow()/ForEachBlock() argument packs
    // into a single capturable object, avoiding NVHPC restrictions on
    // (a) extended lambdas inside generic lambdas and (b) pack element capture.

    template<std::size_t I, typename T>
    struct DTupleElem
    {
      T val;
    };

    template<typename Seq, typename... Ts>
    struct DTupleBase;

    template<std::size_t... Is, typename... Ts>
    struct DTupleBase<std::index_sequence<Is...>, Ts...> : DTupleElem<Is, Ts>...
    {
      KOKKOS_DEFAULTED_FUNCTION DTupleBase() = default;
      template<typename... Us>
      KOKKOS_INLINE_FUNCTION explicit DTupleBase(Us&&... us) : DTupleElem<Is, Ts>{ std::forward<Us>(us) }...
      {
      }
    };

    template<typename... Ts>
    struct DeviceTuple : DTupleBase<std::index_sequence_for<Ts...>, Ts...>
    {
      using Base = DTupleBase<std::index_sequence_for<Ts...>, Ts...>;
      using Base::Base;
      static constexpr std::size_t N = sizeof...(Ts);
    };

    /// @brief TypeAt: get the I-th type from a pack without std::tuple_element
    template<std::size_t I, typename T, typename...>
    struct TypeAtHelper
    {
      using type = T;
    };
    template<std::size_t I, typename T, typename... Rest>
    struct TypeAt : TypeAt<I - 1, Rest...>
    {
    };
    template<typename T, typename... Rest>
    struct TypeAt<0, T, Rest...> : TypeAtHelper<0, T>
    {
    };

    template<std::size_t I, typename... Ts>
    KOKKOS_INLINE_FUNCTION auto& dt_get(DeviceTuple<Ts...>& t) noexcept
    {
      using T = typename TypeAt<I, Ts...>::type;
      return static_cast<DTupleElem<I, T>&>(t).val;
    }

    template<std::size_t I, typename... Ts>
    KOKKOS_INLINE_FUNCTION const auto& dt_get(const DeviceTuple<Ts...>& t) noexcept
    {
      using T = typename TypeAt<I, Ts...>::type;
      return static_cast<const DTupleElem<I, T>&>(t).val;
    }

    template<typename... Ts>
    KOKKOS_INLINE_FUNCTION DeviceTuple<std::decay_t<Ts>...> make_device_tuple(Ts&&... ts)
    {
      return DeviceTuple<std::decay_t<Ts>...>(std::forward<Ts>(ts)...);
    }

  }  // namespace detail

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
    class GroupView;
    class ConstGroupView;
    using ViewType = GroupView;
    using ConstViewType = ConstGroupView;
    using KokkosViewType = Kokkos::View<T*>;
    using HostViewType = typename KokkosViewType::host_mirror_type;
    using TeamPolicyType = Kokkos::TeamPolicy<>;
    using TeamMember = typename TeamPolicyType::member_type;
    template<class VecT>
    using VectorType = KokkosPaddedVector<VecT, L>;

   private:
    /// Device-side (or unified) view — the Kokkos mirror of VectorMatrix::data_
    KokkosViewType view_;

   public:
    // -----------------------------------------------------------------------
    // CUDA-visible implementation types
    // Public so they can appear in __global__ kernel template arguments, which
    // CUDA forbids for private/protected nested types.  These are internal
    // implementation details and are not part of the stable user-facing API.
    // -----------------------------------------------------------------------

    /// @brief Device-safe handle for a mutable KokkosDenseMatrix argument to
    ///        Function()/ForEachRow().
    ///
    /// Only the flat Kokkos::View and column count are captured -- both are trivially
    /// copyable -- so this handle (rather than the matrix object itself, which owns a
    /// non-trivially-copyable host std::vector) is what gets captured by value into a
    /// device lambda.  Only ever used internally by MakeHandle()/BuildGroupView() below.
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

   private:
    /// @brief Build a device-safe handle for one Function()/ForEachRow() argument.
    ///
    /// Kokkos matrix arguments are reduced to their (trivially copyable) View + column
    /// count.
    template<typename Arg>
    static auto MakeHandle(Arg&& arg)
    {
      using ArgType = std::remove_reference_t<Arg>;
      if constexpr (KokkosVectorLike<std::remove_cvref_t<ArgType>>)
      {
        return arg.GetView();
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

    /// @brief Construct the appropriate GroupView/ConstGroupView 
    ///        Runs on-device (called from within a KOKKOS_LAMBDA).
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
        return std::forward<Handle>(handle);
      }
    }

   public:
    /// @brief Kokkos functor for dispatching Function() over complete groups.
    ///        Avoids NVHPC restrictions on extended lambdas inside generic lambdas
    ///        and parameter-pack capture in device lambdas.
    template<typename Func, typename HandlesTuple>
    struct FunctionMainFunctor
    {
      Func func_;
      HandlesTuple handles_;

      template<std::size_t... Is>
      KOKKOS_INLINE_FUNCTION void dispatch(Index group, const TeamMember& team, std::index_sequence<Is...>) const
      {
        func_(BuildGroupView(detail::dt_get<Is>(handles_), group, L, team)...);
      }

      KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team) const
      {
        dispatch(static_cast<Index>(team.league_rank()), team, std::make_index_sequence<HandlesTuple::N>{});
      }
    };

    /// @brief Kokkos functor for dispatching Function() over the tail group.
    template<typename Func, typename HandlesTuple>
    struct FunctionTailFunctor
    {
      Func func_;
      HandlesTuple handles_;
      Index num_complete_groups_;
      Index remaining_;

      template<std::size_t... Is>
      KOKKOS_INLINE_FUNCTION void dispatch(const TeamMember& team, std::index_sequence<Is...>) const
      {
        func_(BuildGroupView(detail::dt_get<Is>(handles_), num_complete_groups_, remaining_, team)...);
      }

      KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team) const
      {
        dispatch(team, std::make_index_sequence<HandlesTuple::N>{});
      }
    };

    /// @brief Kokkos functor for dispatching ForEachRow() via RangePolicy (L == 1).
    template<typename Func, typename ArgsTuple>
    struct ForEachRowRangeFunctor
    {
      Func func_;
      KokkosViewType view_;
      Index y_dim_;
      ArgsTuple args_;

      template<std::size_t... Is>
      KOKKOS_INLINE_FUNCTION void dispatch(Index row, std::index_sequence<Is...>) const
      {
        func_(KokkosDenseMatrix<T, L>::GetTopLevelRowElement(view_, y_dim_, row, detail::dt_get<Is>(args_))...);
      }

      KOKKOS_INLINE_FUNCTION void operator()(Index row) const
      {
        dispatch(row, std::make_index_sequence<ArgsTuple::N>{});
      }
    };

    /// @brief Kokkos functor for dispatching ForEachRow() over complete groups via TeamPolicy.
    template<typename Func, typename ArgsTuple>
    struct ForEachRowTeamFunctor
    {
      Func func_;
      KokkosViewType view_;
      Index y_dim_;
      ArgsTuple args_;

      template<std::size_t... Is>
      KOKKOS_INLINE_FUNCTION void dispatch(const TeamMember& team, Index group, std::index_sequence<Is...>) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, L),
            [&](const Index row_in_group)
            {
              const Index row = group * L + row_in_group;
              func_(KokkosDenseMatrix<T, L>::GetTopLevelRowElement(view_, y_dim_, row, detail::dt_get<Is>(args_))...);
            });
      }

      KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team) const
      {
        dispatch(team, static_cast<Index>(team.league_rank()), std::make_index_sequence<ArgsTuple::N>{});
      }
    };

    /// @brief Kokkos functor for dispatching ForEachRow() over the tail group via TeamPolicy.
    template<typename Func, typename ArgsTuple>
    struct ForEachRowTailFunctor
    {
      Func func_;
      KokkosViewType view_;
      Index y_dim_;
      ArgsTuple args_;
      Index num_complete_groups_;
      Index remaining_;

      template<std::size_t... Is>
      KOKKOS_INLINE_FUNCTION void dispatch(const TeamMember& team, std::index_sequence<Is...>) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, remaining_),
            [&](const Index row_in_group)
            {
              const Index row = num_complete_groups_ * L + row_in_group;
              func_(KokkosDenseMatrix<T, L>::GetTopLevelRowElement(view_, y_dim_, row, detail::dt_get<Is>(args_))...);
            });
      }

      KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team) const
      {
        dispatch(team, std::make_index_sequence<ArgsTuple::N>{});
      }
    };

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
      auto h_view = Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
          this->data_.data(), this->data_.size());
      Kokkos::deep_copy(view_, h_view);
    }

    /// @brief Copy device view data back to host (MICM's data_)
    ///
    /// TODO: Move host view to class member
    void CopyToHost()
    {
      auto h_view = Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
          this->data_.data(), this->data_.size());
      Kokkos::deep_copy(h_view, view_);
    }

    KokkosViewType GetView() const
    {
      return view_;
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

    /// @brief Set every element on the device to a given value
    void Fill(T val)
    {
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
      KokkosViewType y_view = view_;
      KokkosViewType x_view = x.view_;
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
      KokkosViewType y_view = view_;
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
      KokkosViewType y_view = view_;
      Kokkos::parallel_for(
          "KokkosDenseMatrix::Min",
          Kokkos::RangePolicy<>(0, y_view.extent(0)),
          KOKKOS_LAMBDA(const Index i) { y_view(i) = y_view(i) < x ? y_view(i) : x; });
    }

    /// @brief Copy the device data from the other Kokkos dense matrix into this one
    void Copy(const KokkosDenseMatrix& other)
    {
      Kokkos::deep_copy(view_, other.view_);
    }

    /// @brief Swap the device data from the other Kokkos dense matrix into this one.
    void Swap(KokkosDenseMatrix& other)
    {
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
      KokkosViewType y_view = view_;
      KokkosViewType a_view = a.view_;
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
      KokkosViewType y_view = view_;
      KokkosViewType a_view = a.view_;
      KokkosViewType b_view = b.view_;
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

      template<KokkosVectorLike Arg>
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
      template<KokkosVectorLike Vec>
      KOKKOS_INLINE_FUNCTION void Fill(Vec& vec, T value) const
      {
        const Index start = group_ * L;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_rows_in_group_), [&](const Index i) { vec[start + i] = value; });
        team_.team_barrier();
      }

      /// @brief Copy src into vec.
      template<KokkosVectorLike Vec, GroupedDenseMatrixColumnView Src>
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

      /// @brief Apply a reduction to each row in this group, on-device via team
      ///        parallelism. The user's function receives its column-view /
      ///        row-variable arguments plus a trailing reference to a per-thread
      ///        accumulator, and accumulates into it (e.g. `acc += x*x`,
      ///        `acc = std::max(acc, x)`). The micm reducer type (Sum/Max/LOr/LAnd)
      ///        is translated to the matching Kokkos reducer, which handles the
      ///        inter-thread join and writes the final result back to
      ///        `reducer.reference()`.
      template<typename Reducer, typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void Reduce(Reducer reducer, Func&& func, Args&&... args) const
      {
        using KokkosReducer = typename detail::ToKokkosReducer<Reducer>::type;
        using AccT = typename KokkosReducer::value_type;
        // Kokkos::parallel_reduce writes into the reducer's destination scalar
        // per call, overwriting whatever was there. To make repeated Reduce() calls
        // accumulate into the caller's `reducer.reference()` -- matching the host
        // Matrix/VectorMatrix semantics -- reduce into a per-team scratch and then
        // join into the caller's destination from a single team member.
        AccT local = Reducer::identity();
        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team_, L),
            [&](const Index row_in_group, AccT& acc) { func(GetRowElement(row_in_group, args)..., acc); },
            KokkosReducer(local));
        Kokkos::single(Kokkos::PerTeam(team_), [&]() { Reducer::join(reducer.reference(), local); });
        team_.team_barrier();
      }

      /// @brief Same as Reduce but guaranteed to skip padding rows.
      template<typename Reducer, typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void ReduceStrict(Reducer reducer, Func&& func, Args&&... args) const
      {
        using KokkosReducer = typename detail::ToKokkosReducer<Reducer>::type;
        using AccT = typename KokkosReducer::value_type;
        AccT local = Reducer::identity();
        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team_, num_rows_in_group_),
            [&](const Index row_in_group, AccT& acc) { func(GetRowElement(row_in_group, args)..., acc); },
            KokkosReducer(local));
        Kokkos::single(Kokkos::PerTeam(team_), [&]() { Reducer::join(reducer.reference(), local); });
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
      KokkosViewType view_;
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

      template<KokkosVectorLike Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg[group_ * L + row_in_group];
      }

     public:
      KOKKOS_INLINE_FUNCTION
      GroupView(KokkosViewType view, Index group, Index y_dim, Index num_rows_in_group, const TeamMember& team)
          : view_(view),
            group_(group),
            y_dim_(y_dim),
            num_rows_in_group_(num_rows_in_group),
            team_(team)
      {
      }

      KOKKOS_INLINE_FUNCTION operator ConstGroupView() const
      {
        return ConstGroupView(view_, group_, y_dim_, num_rows_in_group_, team_);
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
      template<KokkosVectorLike Src>
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
      template<KokkosVectorLike Vec>
      KOKKOS_INLINE_FUNCTION void Fill(Vec& vec, T value) const
      {
        const Index start = group_ * L;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_rows_in_group_), [&](const Index i) { vec[start + i] = value; });
        team_.team_barrier();
      }

      /// @brief Copy src into vec.
      template<KokkosVectorLike Vec, GroupedDenseMatrixColumnView Src>
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

      /// @brief Apply a reduction to each row in this group, on-device via team
      ///        parallelism. See ConstGroupView::Reduce for details.
      template<typename Reducer, typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void Reduce(Reducer reducer, Func&& func, Args&&... args) const
      {
        using KokkosReducer = typename detail::ToKokkosReducer<Reducer>::type;
        using AccT = typename KokkosReducer::value_type;
        AccT local = Reducer::identity();
        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team_, L),
            [&](const Index row_in_group, AccT& acc) { func(GetRowElement(row_in_group, args)..., acc); },
            KokkosReducer(local));
        Kokkos::single(Kokkos::PerTeam(team_), [&]() { Reducer::join(reducer.reference(), local); });
        team_.team_barrier();
      }

      /// @brief Same as Reduce but guaranteed to skip padding rows.
      template<typename Reducer, typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void ReduceStrict(Reducer reducer, Func&& func, Args&&... args) const
      {
        using KokkosReducer = typename detail::ToKokkosReducer<Reducer>::type;
        using AccT = typename KokkosReducer::value_type;
        AccT local = Reducer::identity();
        Kokkos::parallel_reduce(
            Kokkos::TeamThreadRange(team_, num_rows_in_group_),
            [&](const Index row_in_group, AccT& acc) { func(GetRowElement(row_in_group, args)..., acc); },
            KokkosReducer(local));
        Kokkos::single(Kokkos::PerTeam(team_), [&]() { Reducer::join(reducer.reference(), local); });
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
              if constexpr (KokkosVectorLike<ArgType>)
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

              if constexpr (KokkosVectorLike<ArgType>)
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

        // Reduce each argument to a device-safe handle before entering device code,
        // then bundle the handles into a DeviceTuple and dispatch via named Kokkos functor
        // structs.  This avoids two NVHPC restrictions on extended (host+device) lambdas:
        //   (a) they cannot be defined inside a generic lambda, and
        //   (b) they cannot capture elements of a parameter pack.
        //
        // Always dispatch via Kokkos::TeamPolicy (never a bare RangePolicy), even when
        // L == 1: GroupView/ConstGroupView require a real TeamMember for
        // TeamThreadRange/team_barrier, and TeamMember is not default-constructible.
        auto dev_handles = detail::make_device_tuple(MakeHandle(invoked_args)...);
        using DH = decltype(dev_handles);

        if (num_complete_groups > 0)
        {
          TeamPolicyType policy(static_cast<int>(num_complete_groups), Kokkos::AUTO);
          Kokkos::parallel_for(
              "KokkosDenseMatrix::Function",
              policy,
              FunctionMainFunctor<std::decay_t<decltype(func)>, DH>{ func, dev_handles });
        }
        if (remaining > 0)
        {
          TeamPolicyType tail_policy(1, Kokkos::AUTO);
          Kokkos::parallel_for(
              "KokkosDenseMatrix::Function(tail)",
              tail_policy,
              FunctionTailFunctor<std::decay_t<decltype(func)>, DH>{
                  func, dev_handles, num_complete_groups, remaining });
        }
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
      KokkosViewType view = view_;

      // Bundle args into a DeviceTuple so that the Kokkos functor structs below
      // can capture a single object rather than a parameter pack, satisfying NVHPC.
      auto args_tuple = detail::make_device_tuple(args...);
      using AT = decltype(args_tuple);

      if constexpr (L == 1)
      {
        if (num_rows > 0)
        {
          Kokkos::parallel_for(
              "KokkosDenseMatrix::ForEachRow",
              Kokkos::RangePolicy<>(0, num_rows),
              ForEachRowRangeFunctor<std::decay_t<Func>, AT>{ func, view, y_dim, args_tuple });
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
            ForEachRowTeamFunctor<std::decay_t<Func>, AT>{ func, view, y_dim, args_tuple });
      }
      if (remaining > 0)
      {
        TeamPolicyType tail_policy(1, Kokkos::AUTO);
        Kokkos::parallel_for(
            "KokkosDenseMatrix::ForEachRow(tail)",
            tail_policy,
            ForEachRowTailFunctor<std::decay_t<Func>, AT>{
                func, view, y_dim, args_tuple, num_complete_groups, remaining });
      }
    }

   private:
    /// @brief Get an element reference for a row at the (ungrouped) matrix level.
    ///        Used by the matrix-level ForEachRow() override.
    template<DenseMatrixColumnView Arg>
    KOKKOS_INLINE_FUNCTION static decltype(auto) GetTopLevelRowElement(KokkosViewType, Index, Index row, Arg&& arg)
    {
      return arg.Data()[(row / L * arg.YDim() + arg.ColumnIndex()) * L + row % L];
    }

    template<BlockVariableView Arg>
    KOKKOS_INLINE_FUNCTION static decltype(auto) GetTopLevelRowElement(KokkosViewType, Index, Index row, Arg&& arg)
    {
      return arg.Get()[row % L];
    }

    template<KokkosVectorLike Arg>
    KOKKOS_INLINE_FUNCTION static decltype(auto) GetTopLevelRowElement(KokkosViewType, Index, Index row, Arg&& arg)
    {
      return arg[row];
    }
  };
}  // namespace micm
