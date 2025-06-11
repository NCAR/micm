// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
namespace micm
{
  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolverInPlace(
      const SparseMatrixPolicy& matrix,
      typename SparseMatrixPolicy::value_type initial_value)
      : LinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>(
            matrix,
            initial_value,
            [](const SparseMatrixPolicy& m) -> LuDecompositionPolicy
            { return LuDecompositionPolicy::template Create<SparseMatrixPolicy>(m); })
  {
  }

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolverInPlace(
      const SparseMatrixPolicy& matrix,
      typename SparseMatrixPolicy::value_type initial_value,
      const std::function<LuDecompositionPolicy(const SparseMatrixPolicy&)> create_lu_decomp)
      : nLij_(),
        Lij_yj_(),
        nUij_Uii_(),
        Uij_xj_(),
        lu_decomp_(create_lu_decomp(matrix))
  {
    MICM_PROFILE_FUNCTION();

    auto lu = lu_decomp_.template GetLUMatrix<SparseMatrixPolicy>(matrix, initial_value, matrix.NumberOfBlocks(), true);
    for (std::size_t i = 0; i < lu.NumRows(); ++i)
    {
      std::size_t nLij = 0;
      for (std::size_t j = 0; j < i; ++j)
      {
        if (lu.IsZero(i, j))
          continue;
        Lij_yj_.push_back(std::make_pair(lu.VectorIndex(0, i, j), j));
        ++nLij;
      }
      nLij_.push_back(nLij);
    }
    for (std::size_t i = lu.NumRows() - 1; i != static_cast<std::size_t>(-1); --i)
    {
      std::size_t nUij = 0;
      for (std::size_t j = i + 1; j < lu.NumColumns(); ++j)
      {
        if (lu.IsZero(i, j))
          continue;
        Uij_xj_.push_back(std::make_pair(lu.VectorIndex(0, i, j), j));
        ++nUij;
      }
      // There must always be a non-zero element on the diagonal
      nUij_Uii_.push_back(std::make_pair(nUij, lu.VectorIndex(0, i, i)));
    }
  };

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void LinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>::Factor(SparseMatrixPolicy& matrix) const
  {
    MICM_PROFILE_FUNCTION();

    lu_decomp_.template Decompose<SparseMatrixPolicy>(matrix);
  }

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  template<class MatrixPolicy>
    requires(!VectorizableDense<MatrixPolicy> || !VectorizableSparse<SparseMatrixPolicy>)
  inline void LinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>::Solve(
      MatrixPolicy& x,
      const SparseMatrixPolicy& lu_matrix) const
  {
    MICM_PROFILE_FUNCTION();

    for (std::size_t i_cell = 0; i_cell < x.NumRows(); ++i_cell)
    {
      auto x_cell = x[i_cell];
      const std::size_t grid_offset = i_cell * lu_matrix.FlatBlockSize();
      auto& y_cell = x_cell;  // Alias x for consistency with equations, but to reuse memory

      // Forward Substitution
      {
        auto y_elem = y_cell.begin();
        auto Lij_yj = Lij_yj_.begin();
        for (auto& nLij : nLij_)
        {
          for (std::size_t i = 0; i < nLij; ++i)
          {
            *y_elem -= lu_matrix.AsVector()[grid_offset + (*Lij_yj).first] * y_cell[(*Lij_yj).second];
            ++Lij_yj;
          }
          ++y_elem;
        }
      }

      // Backward Substitution
      {
        auto x_elem = std::next(x_cell.end(), -1);
        auto Uij_xj = Uij_xj_.begin();
        for (auto& nUij_Uii : nUij_Uii_)
        {
          // x_elem starts out as y_elem from the previous loop
          for (std::size_t i = 0; i < nUij_Uii.first; ++i)
          {
            *x_elem -= lu_matrix.AsVector()[grid_offset + (*Uij_xj).first] * x_cell[(*Uij_xj).second];
            ++Uij_xj;
          }

          *(x_elem) /= lu_matrix.AsVector()[grid_offset + nUij_Uii.second];
          // don't iterate before the beginning of the vector
          if (x_elem != x_cell.begin())
          {
            --x_elem;
          }
        }
      }
    }
  }

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  template<class MatrixPolicy>
    requires(VectorizableDense<MatrixPolicy> && VectorizableSparse<SparseMatrixPolicy>)
  inline void LinearSolverInPlace<SparseMatrixPolicy, LuDecompositionPolicy>::Solve(
      MatrixPolicy& x,
      const SparseMatrixPolicy& lu_matrix) const
  {
    MICM_PROFILE_FUNCTION();
    constexpr std::size_t n_cells = MatrixPolicy::GroupVectorSize();
    // Loop over groups of blocks
    for (std::size_t i_group = 0; i_group < x.NumberOfGroups(); ++i_group)
    {
      auto x_group = std::next(x.AsVector().begin(), i_group * x.GroupSize());
      auto LU_group = std::next(lu_matrix.AsVector().begin(), i_group * lu_matrix.GroupSize());
      // Forward Substitution
      {
        auto y_elem = x_group;
        auto Lij_yj = Lij_yj_.begin();
        for (auto& nLij : nLij_)
        {
          for (std::size_t i = 0; i < nLij; ++i)
          {
            const std::size_t Lij_yj_first = (*Lij_yj).first;
            const std::size_t Lij_yj_second_times_n_cells = (*Lij_yj).second * n_cells;
            auto LU_group_it = LU_group + Lij_yj_first;
            auto x_group_it = x_group + Lij_yj_second_times_n_cells;
            auto y_elem_it = y_elem;
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              *(y_elem_it++) -= *(LU_group_it++) * *(x_group_it++);
            ++Lij_yj;
          }
          y_elem += n_cells;
        }
      }

      // Backward Substitution
      {
        auto x_elem = std::next(x_group, x.GroupSize() - n_cells);
        auto Uij_xj = Uij_xj_.begin();
        for (auto& nUij_Uii : nUij_Uii_)
        {
          // x_elem starts out as y_elem from the previous loop
          for (std::size_t i = 0; i < nUij_Uii.first; ++i)
          {
            const std::size_t Uij_xj_first = (*Uij_xj).first;
            const std::size_t Uij_xj_second_times_n_cells = (*Uij_xj).second * n_cells;
            auto LU_group_it = LU_group + Uij_xj_first;
            auto x_group_it = x_group + Uij_xj_second_times_n_cells;
            auto x_elem_it = x_elem;
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              *(x_elem_it++) -= *(LU_group_it++) * *(x_group_it++);
            ++Uij_xj;
          }
          const std::size_t nUij_Uii_second = nUij_Uii.second;
          auto LU_group_it = LU_group + nUij_Uii_second;
          auto x_elem_it = x_elem;
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            *(x_elem_it++) /= *(LU_group_it++);

          // don't iterate before the beginning of the vector
          const std::size_t x_elem_distance = std::distance(x.AsVector().begin(), x_elem);
          x_elem -= std::min(n_cells, x_elem_distance);
        }
      }
    }
  }
}  // namespace micm