// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_padded_vector.hpp>
#include <micm/kokkos/util/kokkos_reducers.hpp>
#include <micm/kokkos/util/kokkos_scalar_view.hpp>
#include <micm/kokkos/util/kokkos_team_policy.hpp>
#include <micm/kokkos/util/kokkos_views.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include <Kokkos_Core.hpp>
#include <cmath>
#include <vector>

namespace micm
{
  /// @brief Provides a Kokkos implementation to the SparseMatrix functionality.
  ///
  /// Inherits from SparseMatrix (the MICM host-side data layout) and maintains
  /// a Kokkos::View as a device-side mirror. The caller must explicitly call
  /// CopyToDevice() / CopyToHost() to synchronize, matching the CUDA matrix pattern.
  template<class T = double, class OrderingPolicy = SparseMatrixVectorOrdering<detail::MICM_KOKKOS_DEFAULT_TEAM_SIZE>>
  class KokkosSparseMatrix : public SparseMatrix<T, OrderingPolicy>
  {
   public:
    static constexpr Index GroupVectorSize()
    {
      return OrderingPolicy::GroupVectorSize();
    }
    using value_type = T;
    template<class U>
    class GroupView;
    using ViewType = GroupView<T>;
    using ConstViewType = GroupView<const T>;
    using KokkosViewType = Kokkos::View<T*>;
    using HostViewType = typename KokkosViewType::host_mirror_type;
    using TeamPolicyType = Kokkos::TeamPolicy<>;
    using TeamMember = typename TeamPolicyType::member_type;
    template<class VecT>
    using VectorType = KokkosPaddedVector<VecT, OrderingPolicy::GroupVectorSize()>;
    template<class ScaT>
    using ScalarType = KokkosScalarView<ScaT>;
    template<class U>
    using SumType = KokkosSum<U>;
    template<class U>
    using MaxType = KokkosMax<U>;
    using LOrType = KokkosLOr;
    using LAndType = KokkosLAnd;

   private:
    /// Number of blocks grouped into a team for on-device iteration.
    static constexpr Index L = OrderingPolicy::GroupVectorSize();

    /// Device-side (or unified) view — the Kokkos mirror of MICM's data_
    KokkosViewType view_;
    /// Host view for copying to and from device
    Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> host_view_;
    Kokkos::View<Index*> diagonal_offsets_;
    bool diagonal_offsets_valid_ = false;

   public:
    // -----------------------------------------------------------------------
    // CUDA-visible implementation types
    // Public so they can appear in __global__ kernel template arguments, which
    // CUDA forbids for private/protected nested types.
    // -----------------------------------------------------------------------

    /// @brief Device-safe handle for a mutable KokkosSparseMatrix argument to
    ///        Function()/ForEachBlock().
    struct SparseMatrixHandle
    {
      Kokkos::View<T*> view_;
      Index flat_block_size_;
    };

    /// @brief Const variant of SparseMatrixHandle. See SparseMatrixHandle for details.
    struct ConstSparseMatrixHandle
    {
      Kokkos::View<const T*> view_;
      Index flat_block_size_;
    };

    /// @brief Device-safe handle for a KokkosDenseMatrix argument mixed into a call to
    ///        Function().
    ///
    /// See SparseMatrixHandle for why only the view + column count are captured.
    struct DenseMatrixArgHandle
    {
      Kokkos::View<T*> view_;
      Index y_dim_;
    };

    /// @brief Const variant of DenseMatrixArgHandle. See DenseMatrixArgHandle for
    ///        details.
    struct ConstDenseMatrixArgHandle
    {
      Kokkos::View<const T*> view_;
      Index y_dim_;
    };

   private:
    /// @brief Build a device-safe handle for one Function()/ForEachBlock() argument.
    template<typename Arg>
    static auto MakeHandle(Arg&& arg)
    {
      using ArgType = std::remove_reference_t<Arg>;
      if constexpr (KokkosVectorLike<std::remove_cvref_t<ArgType>>)
      {
        return arg.GetView();
      }
      else if constexpr (SparseMatrixConcept<std::remove_cvref_t<ArgType>>)
      {
        if constexpr (std::is_const_v<ArgType>)
        {
          return ConstSparseMatrixHandle{ arg.GetView(), arg.FlatBlockSize() };
        }
        else
        {
          return SparseMatrixHandle{ arg.GetView(), arg.FlatBlockSize() };
        }
      }
      else
      {
        // A KokkosDenseMatrix argument (e.g. state variables).
        if constexpr (std::is_const_v<ArgType>)
        {
          return ConstDenseMatrixArgHandle{ arg.GetView(), arg.NumColumns() };
        }
        else
        {
          return DenseMatrixArgHandle{ arg.GetView(), arg.NumColumns() };
        }
      }
    }

    /// @brief Construct the appropriate GroupView/ConstGroupView (or the matching
    ///        KokkosDenseMatrix group view) for one handle produced by MakeHandle().
    ///        Runs on-device (called from within a KOKKOS_LAMBDA).
    template<typename Handle>
    KOKKOS_INLINE_FUNCTION static decltype(auto)
    BuildGroupView(Handle&& handle, Index group, Index count, const TeamMember& team)
    {
      using HandleType = std::remove_cvref_t<Handle>;
      if constexpr (std::is_same_v<HandleType, SparseMatrixHandle>)
      {
        return GroupView<T>(handle.view_, group, handle.flat_block_size_, count, team);
      }
      else if constexpr (std::is_same_v<HandleType, ConstSparseMatrixHandle>)
      {
        return GroupView<const T>(handle.view_, group, handle.flat_block_size_, count, team);
      }
      else if constexpr (std::is_same_v<HandleType, DenseMatrixArgHandle>)
      {
        return typename KokkosDenseMatrix<T, L>::template GroupView<T>(handle.view_, group, handle.y_dim_, count, team);
      }
      else if constexpr (std::is_same_v<HandleType, ConstDenseMatrixArgHandle>)
      {
        return
            typename KokkosDenseMatrix<T, L>::template GroupView<const T>(handle.view_, group, handle.y_dim_, count, team);
      }
      else if (KokkosVectorLike<HandleType>)
      {
        return std::forward<Handle>(handle);
      }
    }

   public:
    /// @brief Kokkos functor for dispatching Function() over complete groups (size L).
    ///        See KokkosDenseMatrix::FunctionMainFunctor for rationale.
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

    /// @brief Kokkos functor for dispatching ForEachBlock() via RangePolicy (L == 1).
    template<typename Func, typename ArgsTuple>
    struct ForEachBlockRangeFunctor
    {
      Func func_;
      KokkosViewType view_;
      Index flat_block_size_;
      ArgsTuple args_;

      template<Index... Is>
      KOKKOS_INLINE_FUNCTION void Dispatch(Index block, std::index_sequence<Is...>) const
      {
        func_(
            KokkosSparseMatrix<T, OrderingPolicy>::GetTopLevelBlockElement(
                view_, flat_block_size_, block, detail::DeviceTupleGet<Is>(args_))...);
      }

      KOKKOS_INLINE_FUNCTION void operator()(Index block) const
      {
        Dispatch(block, std::make_index_sequence<ArgsTuple::N>{});
      }
    };

    /// @brief Kokkos functor for dispatching ForEachBlock() over complete groups via TeamPolicy.
    template<typename Func, typename ArgsTuple>
    struct ForEachBlockTeamFunctor
    {
      Func func_;
      KokkosViewType view_;
      Index flat_block_size_;
      ArgsTuple args_;

      template<Index... Is>
      KOKKOS_INLINE_FUNCTION void Dispatch(const TeamMember& team, Index group, std::index_sequence<Is...>) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, L),
            [&](const Index block_in_group)
            {
              const Index block = group * L + block_in_group;
              func_(
                  KokkosSparseMatrix<T, OrderingPolicy>::GetTopLevelBlockElement(
                      view_, flat_block_size_, block, detail::DeviceTupleGet<Is>(args_))...);
            });
      }

      KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team) const
      {
        Dispatch(team, static_cast<Index>(team.league_rank()), std::make_index_sequence<ArgsTuple::N>{});
      }
    };

    /// @brief Kokkos functor for dispatching ForEachBlock() over the tail group via TeamPolicy.
    template<typename Func, typename ArgsTuple>
    struct ForEachBlockTailFunctor
    {
      Func func_;
      KokkosViewType view_;
      Index flat_block_size_;
      ArgsTuple args_;
      Index num_complete_groups_;
      Index remaining_;

      template<Index... Is>
      KOKKOS_INLINE_FUNCTION void Dispatch(const TeamMember& team, std::index_sequence<Is...>) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, remaining_),
            [&](const Index block_in_group)
            {
              const Index block = num_complete_groups_ * L + block_in_group;
              func_(
                  KokkosSparseMatrix<T, OrderingPolicy>::GetTopLevelBlockElement(
                      view_, flat_block_size_, block, detail::DeviceTupleGet<Is>(args_))...);
            });
      }

      KOKKOS_INLINE_FUNCTION void operator()(const TeamMember& team) const
      {
        Dispatch(team, std::make_index_sequence<ArgsTuple::N>{});
      }
    };

    KokkosSparseMatrix()
        : SparseMatrix<T, OrderingPolicy>()
    {
    }

    KokkosSparseMatrix(const SparseMatrixBuilder<T, OrderingPolicy>& builder, bool indexing_only = false)
        : SparseMatrix<T, OrderingPolicy>(builder, indexing_only),
          view_("sparse_matrix", this->data_.size()),
          host_view_(this->data_.data(), this->data_.size())
    {
    }

    KokkosSparseMatrix& operator=(const SparseMatrixBuilder<T, OrderingPolicy>& builder)
    {
      SparseMatrix<T, OrderingPolicy>::operator=(builder);
      diagonal_offsets_ = Kokkos::View<Index*>();
      diagonal_offsets_valid_ = false;
      view_ = KokkosViewType("sparse_matrix", this->data_.size());
      host_view_ = Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
          this->data_.data(), this->data_.size());
      return *this;
    }

    KokkosSparseMatrix(const KokkosSparseMatrix& other)
        : SparseMatrix<T, OrderingPolicy>(other),
          view_("spase_matrix", other.view_.extent(0)),
          host_view_(this->data_.data(), this->data_.size())
    {
      Kokkos::deep_copy(view_, other.view_);
    }

    KokkosSparseMatrix& operator=(const KokkosSparseMatrix& other)
    {
      if (this == &other)
      {
        return *this;
      }
      SparseMatrix<T, OrderingPolicy>::operator=(other);
      diagonal_offsets_ = Kokkos::View<Index*>();
      diagonal_offsets_valid_ = false;
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
    /// @param init Initial value for scalar
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
          "KokkosSparseMatrix::Fill", Kokkos::RangePolicy<>(0, fill_view.extent(0)), KOKKOS_LAMBDA(const Index i) {
            fill_view(i) = val;
          });
    }

    KokkosSparseMatrix& operator=(T val)
    {
      Fill(val);
      return *this;
    }

    /// @brief Add a value to every diagonal element of every block, on-device.
    void AddToDiagonal(T value)
    {
      if (!diagonal_offsets_valid_)
      {
        const std::vector<Index> offsets = this->DiagonalIndices(0);
        diagonal_offsets_ = Kokkos::View<Index*>("diagonal_offsets", offsets.size());
        if (!offsets.empty())
        {
          auto h_offsets = Kokkos::View<const Index*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
              offsets.data(), offsets.size());
          Kokkos::deep_copy(diagonal_offsets_, h_offsets);
        }
        diagonal_offsets_valid_ = true;
      }
      const auto num_diagonal_elements = static_cast<Index>(diagonal_offsets_.extent(0));
      if (num_diagonal_elements == 0)
      {
        return;
      }

      KokkosViewType data_view = view_;
      Kokkos::View<Index*> offsets_view = diagonal_offsets_;
      const Index group_size = L * this->FlatBlockSize();
      const Index num_groups = (this->NumberOfBlocks() + L - 1) / L;

      Kokkos::parallel_for(
          "KokkosSparseMatrix::AddToDiagonal",
          Kokkos::RangePolicy<>(0, num_groups * num_diagonal_elements * L),
          KOKKOS_LAMBDA(const Index idx) {
            const Index block_in_group = idx % L;
            const Index rest = idx / L;
            const Index offset_index = rest % num_diagonal_elements;
            const Index group = rest / num_diagonal_elements;
            data_view(group * group_size + offsets_view(offset_index) + block_in_group) += value;
          });
    }

    /// @brief Access the non-zero element at a precomputed flat index (from
    ///        VectorIndex(0, row, col)) in every block, for direct on-device
    ///        modification via ForEachBlock().
    KOKKOS_INLINE_FUNCTION KokkosBlockView<T, L> GetBlockView(  // NOLINT(readability-identifier-naming)
        Index vector_index) const
    {
      return KokkosBlockView<T, L>(view_.data(), this->FlatBlockSize(), vector_index);
    }

    /// @brief Const variant of GetBlockView(vector_index). See GetBlockView for
    ///        details.
    KOKKOS_INLINE_FUNCTION KokkosBlockView<const T, L> GetConstBlockView(Index vector_index) const
    {
      return KokkosConstBlockView<T, L>(view_.data(), this->FlatBlockSize(), vector_index);
    }

    /// @brief Apply a function to each block of the matrix.
    /// @tparam Func The lambda/function type
    /// @tparam Args The types of the block view / block variable / vector arguments
    /// @param func The function to apply to each block
    /// @param args Block views, block variables, or vectors
    template<typename Func, typename... Args>
    void ForEachBlock(Func&& func, Args&&... args)
    {
      const Index num_blocks = this->NumberOfBlocks();
      const Index flat_block_size = this->FlatBlockSize();
      KokkosViewType view = view_;

      // Bundle args into a DeviceTuple (see KokkosDenseMatrix::ForEachRow for rationale).
      auto args_tuple = detail::MakeDeviceTuple(args...);
      using AT = decltype(args_tuple);

      if constexpr (L == 1)
      {
        if (num_blocks > 0)
        {
          Kokkos::parallel_for(
              "KokkosSparseMatrix::ForEachBlock",
              Kokkos::RangePolicy<>(0, num_blocks),
              ForEachBlockRangeFunctor<std::decay_t<Func>, AT>{ func, view, flat_block_size, args_tuple });
        }
        return;
      }

      const auto num_complete_groups = static_cast<Index>(num_blocks / L);
      const Index remaining = num_blocks % L;

      if (num_complete_groups > 0)
      {
        const ForEachBlockTeamFunctor<std::decay_t<Func>, AT> team_functor{ func, view, flat_block_size, args_tuple };
        TeamPolicyType policy(static_cast<int>(num_complete_groups), detail::TeamSizeForL<L>(team_functor));
        Kokkos::parallel_for("KokkosSparseMatrix::ForEachBlock", policy, team_functor);
      }
      if (remaining > 0)
      {
        const ForEachBlockTailFunctor<std::decay_t<Func>, AT> team_functor{
          func, view, flat_block_size, args_tuple, num_complete_groups, remaining
        };
        TeamPolicyType tail_policy(1, detail::TeamSizeForL<L>(team_functor));
        Kokkos::parallel_for("KokkosSparseMatrix::ForEachBlock(tail)", tail_policy, team_functor);
      }
    }

    /// @brief GroupView provides a team-parallel mutable view of a single
    ///        block-group of L blocks for iteration on-device.
    template<class U>
    class GroupView
    {
     public:
      using GroupedBlockView = micm::KokkosGroupedBlockView<U>;
      using GroupedConstBlockView = micm::KokkosGroupedBlockView<const U>;

     private:
      Kokkos::View<U*> view_;
      Index group_;
      Index flat_block_size_;
      Index num_blocks_in_group_;
      TeamMember team_;

      template<SparseMatrixBlockView Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg) const
      {
        return arg.Data()[group_ * arg.FlatBlockSize() * L + arg.ElementPosition() + block_in_group];
      }

      template<GroupedSparseMatrixBlockView Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg) const
      {
        return arg.group_base_[arg.block_offset_ + block_in_group];
      }

      template<GroupedDenseMatrixColumnView Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg) const
      {
        return arg.base_[block_in_group];
      }

      template<BlockVariableView Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg) const
      {
        if constexpr (L > 1)
        {
          return arg.Get()[block_in_group];
        }
        else
        {
          return arg.Get();
        }
      }

      template<KokkosVectorLike Arg>
      KOKKOS_INLINE_FUNCTION decltype(auto) GetBlockElement(Index block_in_group, Arg&& arg) const
      {
        if constexpr (L > 1)
        {
          return arg[group_ * L + block_in_group];
        }
        else
        {
          return arg[group_];
        }
      }

     public:
      KOKKOS_INLINE_FUNCTION
      GroupView(Kokkos::View<U*> view, Index group, Index flat_block_size, Index num_blocks_in_group, const TeamMember& team)
          : view_(std::move(view)),
            group_(group),
            flat_block_size_(flat_block_size),
            num_blocks_in_group_(num_blocks_in_group),
            team_(team)
      {
      }

      // NOLINTNEXTLINE(modernize-use-constraints) nvhpc warnings when constraints are used
      template<class V = U, std::enable_if_t<!std::is_const_v<V>, int> = 0>
      KOKKOS_INLINE_FUNCTION operator GroupView<const T>() const
      {
        return GroupView<const T>(view_, group_, flat_block_size_, num_blocks_in_group_, team_);
      }

      KOKKOS_INLINE_FUNCTION GroupedConstBlockView GetConstBlockView(Index vector_index) const
      {
        return { view_.data() + group_ * flat_block_size_ * L, vector_index };
      }

      KOKKOS_INLINE_FUNCTION GroupedBlockView GetBlockView(Index vector_index) const
      {
        return { view_.data() + group_ * flat_block_size_ * L, vector_index };
      }

      KOKKOS_INLINE_FUNCTION KokkosBlockVariable<T, L> GetBlockVariable() const
      {
        return KokkosBlockVariable<T, L>();
      }

      /// @brief Assign value to every cell of a grouped block view.
      KOKKOS_INLINE_FUNCTION void Fill(GroupedBlockView view, const T value) const
      {
        T* dst = view.group_base_ + view.block_offset_;
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_, L), [&](const Index i) { dst[i] = value; });
        team_.team_barrier();
      }

      /// @brief Copy src into dst_view.
      template<GroupedSparseMatrixBlockView Src>
      KOKKOS_INLINE_FUNCTION void Copy(GroupedBlockView dst_view, Src&& src_view) const
      {
        T* dst = dst_view.group_base_ + dst_view.block_offset_;
        const T* src = src_view.group_base_ + src_view.block_offset_;
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_, L), [&](const Index i) { dst[i] = src[i]; });
        team_.team_barrier();
      }

      /// @brief Copy src into dst_view.
      template<KokkosVectorLike Src>
      KOKKOS_INLINE_FUNCTION void Copy(GroupedBlockView dst_view, Src&& src) const
      {
        T* dst = dst_view.group_base_ + dst_view.block_offset_;
        const Index start = group_ * L;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_blocks_in_group_), [&](const Index i) { dst[i] = src[start + i]; });
        team_.team_barrier();
      }

      /// @brief Assign value to every cell of dst.
      template<BlockVariableView Dst>
      KOKKOS_INLINE_FUNCTION void Fill(Dst&& dst, const T value) const
      {
        auto& storage = dst.Get();
        if constexpr (L > 1)
        {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team_, L), [&](const Index i) { storage[i] = value; });
          team_.team_barrier();
        }
        else
        {
          storage = value;
        }
      }

      /// @brief Copy a sparse-block value into the caller-owned block-variable temp.
      template<BlockVariableView Dst, GroupedSparseMatrixBlockView Src>
      KOKKOS_INLINE_FUNCTION void Copy(Dst&& dst, Src&& src) const
      {
        auto& storage = dst.Get();
        if constexpr (L > 1)
        {
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team_, L),
              [&](const Index i) { storage[i] = src.group_base_[src.block_offset_ + i]; });
          team_.team_barrier();
        }
        else
        {
          storage = src.group_base_[src.block_offset_];
        }
      }

      /// @brief Assign value to every element of vec.
      template<KokkosVectorLike Vec>
      KOKKOS_INLINE_FUNCTION void Fill(Vec& vec, const T value) const
      {
        const Index start = group_ * L;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_blocks_in_group_), [&](const Index i) { vec[start + i] = value; });
        team_.team_barrier();
      }

      /// @brief Copy src into vec.
      template<KokkosVectorLike Vec, GroupedSparseMatrixBlockView Src>
      KOKKOS_INLINE_FUNCTION void Copy(Vec& vec, Src&& src) const
      {
        const Index start = group_ * L;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_blocks_in_group_),
            [&](const Index i) { vec[start + i] = src.group_base_[src.block_offset_ + i]; });
        team_.team_barrier();
      }

      /// @brief Apply the provided function to every block in this group, including
      ///        any trailing padding blocks.
      template<typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void ForEachBlock(Func&& func, Args&&... args) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, L),
            [&](const Index block_in_group) { func(GetBlockElement(block_in_group, args)...); });
        team_.team_barrier();
      }

      /// @brief Same as ForEachBlock but guaranteed to skip padding blocks.
      template<typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void ForEachBlockStrict(Func&& func, Args&&... args) const
      {
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_, num_blocks_in_group_),
            [&](const Index block_in_group) { func(GetBlockElement(block_in_group, args)...); });
        team_.team_barrier();
      }

      /// @brief Apply a reduction to each row in this group, on-device via team
      ///        parallelism. The user's function receives its column-view /
      ///        row-variable arguments plus a trailing reference to a per-thread
      ///        accumulator, and accumulates into it (e.g. `acc += x*x`,
      ///        `acc = std::max(acc, x)`). The micm reducer type (Sum/Max/LOr/LAnd)
      ///        is translated to the matching Kokkos reducer, which handles the
      ///        inter-thread join and writes the final result back to
      ///        `reducer.Reference()`.
      template<typename Reducer, typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void Reduce(Reducer reducer, Func&& func, Args&&... args) const
      {
        using AccT = decltype(Reducer::Identity());
        reducer.TeamReduce(
            team_,
            L,
            [&](const Index block_in_group, AccT& acc)
            { func(GetBlockElement(block_in_group, std::forward<Args>(args))..., acc); });
        team_.team_barrier();
      }

      /// @brief Same as Reduce but guaranteed to skip padding rows.
      template<typename Reducer, typename Func, typename... Args>
      KOKKOS_INLINE_FUNCTION void ReduceStrict(Reducer reducer, Func&& func, Args&&... args) const
      {
        using AccT = decltype(Reducer::Identity());
        reducer.TeamReduce(
            team_,
            num_blocks_in_group_,
            [&](const Index block_in_group, AccT& acc)
            { func(GetBlockElement(block_in_group, std::forward<Args>(args))..., acc); });
        team_.team_barrier();
      }
    };

    /// @brief Create a function that can be applied to Kokkos sparse and dense
    ///        matrices and vectors, executing on-device using team parallelism.
    template<typename Func, typename... Args>
    static auto Function(Func&& func, Args&... args)
    {
      constexpr Index expected_L = OrderingPolicy::GroupVectorSize();
      Index index = 0;
      (
          [&](auto& arg)
          {
            using ArgType = std::remove_cvref_t<decltype(arg)>;
            if constexpr (!KokkosVectorLike<ArgType>)
            {
              constexpr Index arg_L = GROUP_VECTOR_SIZE_V<ArgType>;
              if (arg_L != expected_L)
              {
                throw MicmException(
                    MICM_ERROR_CATEGORY_MATRIX,
                    MICM_MATRIX_ERROR_CODE_INVALID_VECTOR,
                    "Incompatible matrix orderings: Matrix " + std::to_string(index) +
                        " has GroupVectorSize=" + std::to_string(arg_L) + " but expected " + std::to_string(expected_L));
              }
            }
            ++index;
          }(args),
          ...);

      auto result = [func = std::forward<Func>(func)](auto&... invoked_args) mutable
      {
        Index num_blocks = 0;
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
                  num_blocks = arg.size();
                  found_first = true;
                }
                else if (arg.size() != num_blocks)
                {
                  throw MicmException(
                      MICM_ERROR_CATEGORY_MATRIX,
                      MICM_MATRIX_ERROR_CODE_INVALID_VECTOR,
                      "Vector size must match block count. Expected " + std::to_string(num_blocks) + " elements but got " +
                          std::to_string(arg.size()));
                }
              }
              else
              {
                Index arg_blocks;
                if constexpr (SparseMatrixConcept<ArgType>)
                {
                  arg_blocks = arg.NumberOfBlocks();
                }
                else
                {
                  arg_blocks = arg.NumRows();
                }

                if (!found_first)
                {
                  num_blocks = arg_blocks;
                  found_first = true;
                }
                else if (arg_blocks != num_blocks)
                {
                  throw MicmException(
                      MICM_ERROR_CATEGORY_MATRIX,
                      MICM_MATRIX_ERROR_CODE_INVALID_VECTOR,
                      "All matrices must have the same number of blocks/rows when invoking function. Expected " +
                          std::to_string(num_blocks) + " but got " + std::to_string(arg_blocks));
                }
              }
              ++idx;
            }(invoked_args),
            ...);

        auto num_complete_groups = static_cast<Index>(num_blocks / L);
        Index remaining = num_blocks % L;

        // Bundle handles into a DeviceTuple and dispatch via named Kokkos functor structs.
        // See KokkosDenseMatrix::Function() for the full rationale on why
        // generic-lambda + KOKKOS_LAMBDA nesting and pack capture must be avoided.
        auto dev_handles = detail::MakeDeviceTuple(MakeHandle(invoked_args)...);
        using DH = decltype(dev_handles);

        if (num_complete_groups > 0)
        {
          const FunctionMainFunctor<std::decay_t<decltype(func)>, DH> team_functor{ func, dev_handles };
          TeamPolicyType policy(static_cast<int>(num_complete_groups), detail::TeamSizeForL<L>(team_functor));
          Kokkos::parallel_for("KokkosSparseMatrix::Function", policy, team_functor);
        }
        if (remaining > 0)
        {
          const FunctionTailFunctor<std::decay_t<decltype(func)>, DH> team_functor{
            func, dev_handles, num_complete_groups, remaining
          };
          TeamPolicyType tail_policy(1, detail::TeamSizeForL<L>(team_functor));
          Kokkos::parallel_for("KokkosSparseMatrix::Function(tail)", tail_policy, team_functor);
        }
      };
      return result;
    }

    /// @brief Get an element reference for a block at the (ungrouped) matrix level.
    ///        Used by the matrix-level ForEachBlock() override.
    template<SparseMatrixBlockView Arg>
    KOKKOS_INLINE_FUNCTION static decltype(auto) GetTopLevelBlockElement(KokkosViewType, Index, Index block, Arg&& arg)
    {
      return arg.Data()[(block / L) * arg.FlatBlockSize() * L + arg.ElementPosition() + block % L];
    }

    template<BlockVariableView Arg>
    KOKKOS_INLINE_FUNCTION static decltype(auto) GetTopLevelBlockElement(KokkosViewType, Index, Index block, Arg&& arg)
    {
      return arg.Get();
    }

    template<KokkosVectorLike Arg>
    KOKKOS_INLINE_FUNCTION static decltype(auto) GetTopLevelBlockElement(KokkosViewType, Index, Index block, Arg&& arg)
    {
      return arg[block];
    }
  };
}  // namespace micm
