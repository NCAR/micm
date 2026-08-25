// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/kokkos/util/kokkos_padded_vector.hpp>
#include <micm/kokkos/util/kokkos_reducers.hpp>
#include <micm/kokkos/util/kokkos_scalar_view.hpp>
#include <micm/kokkos/util/kokkos_team_policy.hpp>
#include <micm/kokkos/util/kokkos_view_category.hpp>
#include <micm/kokkos/util/kokkos_views.hpp>
#include <micm/util/reducers.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <Kokkos_Core.hpp>
#include <vector>

#ifndef MICM_KOKKOS_DEFAULT_VECTOR_SIZE
  #if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
    #define MICM_KOKKOS_DEFAULT_VECTOR_SIZE 32
  #else
    #define MICM_KOKKOS_DEFAULT_VECTOR_SIZE MICM_DEFAULT_VECTOR_SIZE
  #endif
#endif

namespace micm
{
  namespace detail
  {
    // ----------------------------------------------------------------
    // Device-compatible tuple for KOKKOS_LAMBDA-safe variadic dispatch
    // ----------------------------------------------------------------
    // No STL runtime dependencies, works inside KOKKOS_FUNCTION contexts.
    // Used to bundle Function()/ForEachRow()/ForEachBlock() argument packs
    // into a single capturable object, avoiding NVHPC restrictions on
    // (a) extended lambdas inside generic lambdas and (b) pack element capture.

    template<Index I, typename T>
    struct DTupleElem
    {
      T val_;
    };

    template<typename Seq, typename... Ts>
    struct DTupleBase;

    template<Index... Is, typename... Ts>
    struct DTupleBase<std::index_sequence<Is...>, Ts...> : DTupleElem<Is, Ts>...
    {
      KOKKOS_DEFAULTED_FUNCTION DTupleBase() = default;
      template<typename... Us>
      KOKKOS_INLINE_FUNCTION explicit DTupleBase(Us&&... us)
          : DTupleElem<Is, Ts>{ std::forward<Us>(us) }...
      {
      }
    };

    template<typename... Ts>
    struct DeviceTuple : DTupleBase<std::index_sequence_for<Ts...>, Ts...>
    {
      using Base = DTupleBase<std::index_sequence_for<Ts...>, Ts...>;
      using Base::Base;
      static constexpr Index N = sizeof...(Ts);
    };

    template<Index I, typename T>
    KOKKOS_INLINE_FUNCTION auto& DeviceTupleGet(DTupleElem<I, T>& elem) noexcept
    {
      return elem.val_;
    }

    template<Index I, typename T>
    KOKKOS_INLINE_FUNCTION const T& DeviceTupleGet(const DTupleElem<I, T>& elem) noexcept
    {
      return elem.val_;
    }

    template<typename... Ts>
    KOKKOS_INLINE_FUNCTION DeviceTuple<std::decay_t<Ts>...> MakeDeviceTuple(Ts&&... ts)
    {
      return DeviceTuple<std::decay_t<Ts>...>(std::forward<Ts>(ts)...);
    }

    constexpr Index MICM_KOKKOS_DEFAULT_TEAM_SIZE = 128;

  }  // namespace detail

  /// @brief Provides a Kokkos implementation to the VectorMatrix functionality.
  ///
  /// Inherits from VectorMatrix (the MICM host-side data layout) and maintains
  /// a Kokkos::View as a device-side mirror. The caller must explicitly call
  /// CopyToDevice() / CopyToHost() to synchronize, matching the CUDA matrix pattern.
  template<class T, Index L = detail::MICM_KOKKOS_DEFAULT_TEAM_SIZE>
  class KokkosDenseMatrix : public VectorMatrix<T, L>
  {
   public:
    static constexpr Index GroupVectorSize()
    {
      return L;
    }
    using value_type = T;
    template<class U>
    class GroupView;
    using ViewType = GroupView<T>;
    using ConstViewType = GroupView<const T>;
    using HostGroupView = typename VectorMatrix<T, L>::GroupView;
    using ConstHostGroupView = typename VectorMatrix<T, L>::ConstGroupView;
    using KokkosViewType = Kokkos::View<T*>;
    using HostViewType = typename KokkosViewType::host_mirror_type;
    using TeamPolicyType = Kokkos::TeamPolicy<>;
    using TeamMember = typename TeamPolicyType::member_type;
    template<class VecT>
    using VectorType = KokkosPaddedVector<VecT, L>;
    template<class ScaT>
    using ScalarType = KokkosScalarView<ScaT>;
    template<class U>
    using SumType = KokkosSum<U>;
    template<class U>
    using MaxType = KokkosMax<U>;
    using LOrType = KokkosLOr;
    using LAndType = KokkosLAnd;

   private:
    /// Device-side (or unified) view — the Kokkos mirror of VectorMatrix::data_
    KokkosViewType view_;
    /// Host view for copying to and from device
    Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> host_view_;

   public:
    // -----------------------------------------------------------------------
    // CUDA-visible implementation types
    // Public so they can appear in __global__ kernel template arguments, which
    // CUDA forbids for private/protected nested types.
    // -----------------------------------------------------------------------

    /// @brief Device-safe handle for a mutable KokkosDenseMatrix argument to
    ///        Function()/ForEachRow().
    struct DenseMatrixHandle
    {
      Kokkos::View<T*> view_;
      Index y_dim_;
    };

    /// @brief Const variant of DenseMatrixHandle. See DenseMatrixHandle for details.
    struct ConstDenseMatrixHandle
    {
      Kokkos::View<const T*> view_;
      Index y_dim_;
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
    KOKKOS_INLINE_FUNCTION static decltype(auto)
    BuildGroupView(Handle&& handle, Index group, Index count, const TeamMember& team)
    {
      using HandleType = std::remove_cvref_t<Handle>;
      if constexpr (std::is_same_v<HandleType, DenseMatrixHandle>)
      {
        return GroupView<T>(handle.view_, group, handle.y_dim_, count, team);
      }
      else if constexpr (std::is_same_v<HandleType, ConstDenseMatrixHandle>)
      {
        return GroupView<const T>(handle.view_, group, handle.y_dim_, count, team);
      }
      else
      {
        return std::forward<Handle>(handle);
      }
    }

   public:
    /// @brief Kokkos functor for dispatching Function() over complete groups (size L).
    ///        Avoids NVHPC restrictions on extended lambdas inside generic lambdas
    ///        and parameter-pack capture in device lambdas.
    template<typename Func, typename HandlesTuple>
    struct FunctionMainFunctor
    {
      Func func_;
      HandlesTuple handles_;

      template<Index... Is>
      KOKKOS_INLINE_FUNCTION void Dispatch(Index group, const TeamMember& team, std::index_sequence<Is...>) const
      {
        func_(BuildGroupView(detail::DeviceTupleGet<Is>(handles_), group, L, team)...);
      }

      KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team) const
      {
        Dispatch(static_cast<Index>(team.league_rank()), team, std::make_index_sequence<HandlesTuple::N>{});
      }
    };

    /// @brief Kokkos functor for dispatching Function() over the tail group (size < L).
    template<typename Func, typename HandlesTuple>
    struct FunctionTailFunctor
    {
      Func func_;
      HandlesTuple handles_;
      Index num_complete_groups_;
      Index remaining_;

      template<Index... Is>
      KOKKOS_INLINE_FUNCTION void Dispatch(const TeamMember& team, std::index_sequence<Is...>) const
      {
        func_(BuildGroupView(detail::DeviceTupleGet<Is>(handles_), num_complete_groups_, remaining_, team)...);
      }

      KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team) const
      {
        Dispatch(team, std::make_index_sequence<HandlesTuple::N>{});
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

      template<Index... Is>
      KOKKOS_INLINE_FUNCTION void Dispatch(Index row, std::index_sequence<Is...>) const
      {
        func_(KokkosDenseMatrix<T, L>::GetTopLevelRowElement(view_, y_dim_, row, detail::DeviceTupleGet<Is>(args_))...);
      }

      KOKKOS_INLINE_FUNCTION void operator()(Index row) const
      {
        Dispatch(row, std::make_index_sequence<ArgsTuple::N>{});
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

      template<Index... Is>
      KOKKOS_INLINE_FUNCTION void Dispatch(const TeamMember& team, Index group, std::index_sequence<Is...>) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, L),
            [&](const Index row_in_group)
            {
              const Index row = group * L + row_in_group;
              func_(
                  KokkosDenseMatrix<T, L>::GetTopLevelRowElement(view_, y_dim_, row, detail::DeviceTupleGet<Is>(args_))...);
            });
      }

      KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team) const
      {
        Dispatch(team, static_cast<Index>(team.league_rank()), std::make_index_sequence<ArgsTuple::N>{});
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

      template<Index... Is>
      KOKKOS_INLINE_FUNCTION void Dispatch(const TeamMember& team, std::index_sequence<Is...>) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, remaining_),
            [&](const Index row_in_group)
            {
              const Index row = num_complete_groups_ * L + row_in_group;
              func_(
                  KokkosDenseMatrix<T, L>::GetTopLevelRowElement(view_, y_dim_, row, detail::DeviceTupleGet<Is>(args_))...);
            });
      }

      KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team) const
      {
        Dispatch(team, std::make_index_sequence<ArgsTuple::N>{});
      }
    };

    /// @brief Kokkos functor for dispatching 2-arg flat ForEach() over the main range.
    template<typename Func>
    struct ForEachFunctor2
    {
      Func func_;
      KokkosViewType y_view_;
      KokkosViewType a_view_;
      KOKKOS_INLINE_FUNCTION void operator()(const Index i) const
      {
        func_(y_view_(i), a_view_(i));
      }
    };

    /// @brief Kokkos functor for dispatching 2-arg flat ForEach() over the tail range.
    template<typename Func>
    struct ForEachTailFunctor2
    {
      Func func_;
      KokkosViewType y_view_;
      KokkosViewType a_view_;
      Index n_;
      Index l_;
      KOKKOS_INLINE_FUNCTION void operator()(const Index idx) const
      {
        const Index flat = n_ + (idx / l_) * L + (idx % l_);
        func_(y_view_(flat), a_view_(flat));
      }
    };

    /// @brief Kokkos functor for dispatching 3-arg flat ForEach() over the main range.
    template<typename Func>
    struct ForEachFunctor3
    {
      Func func_;
      KokkosViewType y_view_;
      KokkosViewType a_view_;
      KokkosViewType b_view_;
      KOKKOS_INLINE_FUNCTION void operator()(const Index i) const
      {
        func_(y_view_(i), a_view_(i), b_view_(i));
      }
    };

    /// @brief Kokkos functor for dispatching 3-arg flat ForEach() over the tail range.
    template<typename Func>
    struct ForEachTailFunctor3
    {
      Func func_;
      KokkosViewType y_view_;
      KokkosViewType a_view_;
      KokkosViewType b_view_;
      Index n_;
      Index l_;
      KOKKOS_INLINE_FUNCTION void operator()(const Index idx) const
      {
        const Index flat = n_ + (idx / l_) * L + (idx % l_);
        func_(y_view_(flat), a_view_(flat), b_view_(flat));
      }
    };

    KokkosDenseMatrix()
        : VectorMatrix<T, L>()
    {
    }

    KokkosDenseMatrix(Index x_dim, Index y_dim)
        : VectorMatrix<T, L>(x_dim, y_dim),
          view_("dense_matrix", this->data_.size()),
          host_view_(this->data_.data(), this->data_.size())
    {
    }

    KokkosDenseMatrix(Index x_dim, Index y_dim, T initial_value)
        : VectorMatrix<T, L>(x_dim, y_dim, initial_value),
          view_("dense_matrix", this->data_.size()),
          host_view_(this->data_.data(), this->data_.size())
    {
      Kokkos::deep_copy(view_, initial_value);
    }

    KokkosDenseMatrix(const std::vector<std::vector<T>>& other)
        : VectorMatrix<T, L>(other),
          view_("dense_matrix", this->data_.size()),
          host_view_(this->data_.data(), this->data_.size())
    {
      CopyToDevice();
    }

    KokkosDenseMatrix(const KokkosDenseMatrix& other)
        : VectorMatrix<T, L>(other),
          view_("dense_matrix", other.view_.extent(0)),
          host_view_(this->data_.data(), this->data_.size())
    {
      Kokkos::deep_copy(view_, other.view_);
    }

    KokkosDenseMatrix& operator=(const KokkosDenseMatrix& other)
    {
      if (this == &other)
      {
        return *this;
      }
      VectorMatrix<T, L>::operator=(other);
      Kokkos::realloc(view_, other.view_.extent(0));
      Kokkos::deep_copy(view_, other.view_);
      host_view_ = Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
          this->data_.data(), this->data_.size());
      return *this;
    }

    /// @brief Copy host data (MICM's data_) to the device view
    void CopyToDevice()
    {
      Kokkos::deep_copy(view_, host_view_);
    }

    /// @brief Copy device view data back to host (MICM's data_)
    void CopyToHost()
    {
      Kokkos::deep_copy(host_view_, view_);
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

    /// @brief Creates a scalar usable with this matrix type in Function lambda captures
    /// @param init initial value for scalar
    /// @return scalar usable in Function() lambda captures
    template<class ScaT>
    ScalarType<ScaT> CompatibleScalar(ScaT init = ScaT{}) const
    {
      return ScalarType<ScaT>(init);
    }

    /// @brief Set every element on the device to a given value
    void Fill(T val)
    {
      KokkosViewType fill_view = view_;
      Kokkos::parallel_for(
          "KokkosDenseMatrix::Fill", Kokkos::RangePolicy<>(0, fill_view.extent(0)), KOKKOS_LAMBDA(const Index i) {
            fill_view(i) = val;
          });
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
      const Index n = (this->NumRows() / L) * L * y_dim;
      Kokkos::parallel_for(
          "KokkosDenseMatrix::Axpy", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const Index i) {
            y_view(i) += alpha * x_view(i);
          });
      const Index l = this->NumRows() % L;
      if (l > 0)
      {
        Kokkos::parallel_for(
            "KokkosDenseMatrix::Axpy(tail)", Kokkos::RangePolicy<>(0, y_dim * l), KOKKOS_LAMBDA(const Index idx) {
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
          "KokkosDenseMatrix::Max", Kokkos::RangePolicy<>(0, y_view.extent(0)), KOKKOS_LAMBDA(const Index i) {
            y_view(i) = y_view(i) > x ? y_view(i) : x;
          });
    }

    /// @brief For each element of the matrix, perform y = min(y, x), where x is a scalar constant.
    ///
    /// Touches every stored cell, including any trailing padding cells.
    void Min(const T& x)
    {
      KokkosViewType y_view = view_;
      Kokkos::parallel_for(
          "KokkosDenseMatrix::Min", Kokkos::RangePolicy<>(0, y_view.extent(0)), KOKKOS_LAMBDA(const Index i) {
            y_view(i) = y_view(i) < x ? y_view(i) : x;
          });
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
      std::swap(host_view_, other.host_view_);
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
      const Index n = static_cast<Index>(this->NumRows() / L) * L * y_dim;
      Kokkos::parallel_for(
          "KokkosDenseMatrix::ForEach",
          Kokkos::RangePolicy<>(0, n),
          ForEachFunctor2<std::decay_t<Func>>{ f, y_view, a_view });
      const Index l = this->NumRows() % L;
      if (l > 0)
      {
        Kokkos::parallel_for(
            "KokkosDenseMatrix::ForEach(tail)",
            Kokkos::RangePolicy<>(0, y_dim * l),
            ForEachTailFunctor2<std::decay_t<Func>>{ f, y_view, a_view, n, l });
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
      const Index n = static_cast<Index>(this->NumRows() / L) * L * y_dim;
      Kokkos::parallel_for(
          "KokkosDenseMatrix::ForEach",
          Kokkos::RangePolicy<>(0, n),
          ForEachFunctor3<std::decay_t<Func>>{ f, y_view, a_view, b_view });
      const Index l = this->NumRows() % L;
      if (l > 0)
      {
        Kokkos::parallel_for(
            "KokkosDenseMatrix::ForEach(tail)",
            Kokkos::RangePolicy<>(0, y_dim * l),
            ForEachTailFunctor3<std::decay_t<Func>>{ f, y_view, a_view, b_view, n, l });
      }
    }

    /// @brief GroupView provides a team-parallel view of a single group of L rows for
    ///        iteration on-device.
    template<class U>
    class GroupView
    {
     public:
      using GroupedColumnView = micm::KokkosGroupedColumnView<U>;
      using GroupedConstColumnView = micm::KokkosGroupedColumnView<const U>;

     private:
      Kokkos::View<U*> view_;
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
        auto& storage = arg.Get();
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(storage)>, T>)
        {
          return storage;  // KokkosBlockVariable<T,1>: scalar
        }
        else
        {
          return storage[row_in_group];  // KokkosRowVariable (always array) or L>1
        }
      }

      template<KokkosVectorLike Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetRowElement(Index row_in_group, Arg&& arg) const
      {
        return arg[group_ * L + row_in_group];
      }

     public:
      KOKKOS_INLINE_FUNCTION
      GroupView(Kokkos::View<U*> view, Index group, Index y_dim, Index num_rows_in_group, const TeamMember& team)
          : view_(std::move(view)),
            group_(group),
            y_dim_(y_dim),
            num_rows_in_group_(num_rows_in_group),
            team_(team)
      {
      }

      // NOLINTNEXTLINE(modernize-use-constraints) nvhpc warnings when constraints are used
      template<class V = U, std::enable_if_t<!std::is_const_v<V>, int> = 0>
      KOKKOS_INLINE_FUNCTION operator GroupView<const T>() const
      {
        return GroupView<const T>(view_, group_, y_dim_, num_rows_in_group_, team_);
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
      KOKKOS_INLINE_FUNCTION void Fill(GroupedColumnView view, const T value) const
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
      KOKKOS_INLINE_FUNCTION void Fill(Dst&& dst, const T value) const
      {
        auto& storage = dst.Get();
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(storage)>, T>)
        {
          storage = value;  // KokkosBlockVariable<T,1>: scalar
        }
        else if constexpr (L == 1)
        {
          storage[0] = value;  // KokkosRowVariable<T,1>: array of 1
        }
        else
        {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team_, L), [&](const Index i) { storage[i] = value; });
        }
        team_.team_barrier();
      }

      /// @brief Copy src into dst.
      template<BlockVariableView Dst, GroupedDenseMatrixColumnView Src>
      KOKKOS_INLINE_FUNCTION void Copy(Dst&& dst, Src&& src) const
      {
        auto& storage = dst.Get();
        if constexpr (std::is_same_v<std::remove_reference_t<decltype(storage)>, T>)
        {
          storage = src.base_[0];  // KokkosBlockVariable<T,1>: scalar
        }
        else if constexpr (L == 1)
        {
          storage[0] = src.base_[0];  // KokkosRowVariable<T,1>: array of 1
        }
        else
        {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team_, L), [&](const Index i) { storage[i] = src.base_[i]; });
        }
        team_.team_barrier();
      }

      /// @brief Assign value to all vec elements.
      template<KokkosVectorLike Vec>
      KOKKOS_INLINE_FUNCTION void Fill(Vec& vec, const T value) const
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
      KOKKOS_INLINE_FUNCTION void Reduce(const Reducer& reducer, Func&& func, Args&&... args) const
      {
        using AccT = decltype(Reducer::Identity());
        reducer.TeamReduce(
            team_,
            L,
            [&](const Index row_in_group, AccT& acc)
            { func(GetRowElement(row_in_group, std::forward<Args>(args))..., acc); });
        team_.team_barrier();
      }

      /// @brief Same as Reduce but guaranteed to skip padding rows.
      template<typename Reducer, typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void ReduceStrict(const Reducer& reducer, Func&& func, Args&&... args) const
      {
        using AccT = decltype(Reducer::Identity());
        reducer.TeamReduce(
            team_,
            num_rows_in_group_,
            [&](const Index row_in_group, AccT& acc)
            { func(GetRowElement(row_in_group, std::forward<Args>(args))..., acc); });
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
        Bool found_first = false;
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

        auto num_complete_groups = static_cast<Index>(num_rows / L);
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
        auto dev_handles = detail::MakeDeviceTuple(MakeHandle(invoked_args)...);
        using DH = decltype(dev_handles);

        if (num_complete_groups > 0)
        {
          const FunctionMainFunctor<std::decay_t<decltype(func)>, DH> team_functor{ func, dev_handles };
          TeamPolicyType policy(static_cast<int>(num_complete_groups), detail::TeamSizeForL<L>(team_functor));
          Kokkos::parallel_for("KokkosDenseMatrix::Function", policy, team_functor);
        }
        if (remaining > 0)
        {
          const FunctionTailFunctor<std::decay_t<decltype(func)>, DH> team_functor{
            func, dev_handles, num_complete_groups, remaining
          };
          TeamPolicyType tail_policy(1, detail::TeamSizeForL<L>(team_functor));
          Kokkos::parallel_for("KokkosDenseMatrix::Function(tail)", tail_policy, team_functor);
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
      auto args_tuple = detail::MakeDeviceTuple(args...);
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

      const auto num_complete_groups = num_rows / L;
      const Index remaining = num_rows % L;

      if (num_complete_groups > 0)
      {
        const ForEachRowTeamFunctor<std::decay_t<Func>, AT> team_functor{ func, view, y_dim, args_tuple };
        TeamPolicyType policy(static_cast<int>(num_complete_groups), detail::TeamSizeForL<L>(team_functor));
        Kokkos::parallel_for("KokkosDenseMatrix::ForEachRow", policy, team_functor);
      }
      if (remaining > 0)
      {
        const ForEachRowTailFunctor<std::decay_t<Func>, AT> team_functor{
          func, view, y_dim, args_tuple, num_complete_groups, remaining
        };
        TeamPolicyType tail_policy(1, detail::TeamSizeForL<L>(team_functor));
        Kokkos::parallel_for("KokkosDenseMatrix::ForEachRow(tail)", tail_policy, team_functor);
      }
    }

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
