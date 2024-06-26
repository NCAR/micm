// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{
  template<class MatrixPolicy>
  inline std::vector<std::size_t> DiagonalMarkowitzReorder(const MatrixPolicy& matrix)
  {
    MICM_PROFILE_FUNCTION();
    const std::size_t order = matrix.NumRows();
    std::vector<std::size_t> perm(order);
    for (std::size_t i = 0; i < order; ++i)
      perm[i] = i;
    MatrixPolicy pattern = matrix;
    for (std::size_t row = 0; row < (order - 1); ++row)
    {
      std::size_t beta = std::pow((order - 1), 2);
      std::size_t max_row = row;
      for (std::size_t col = row; col < order; ++col)
      {
        std::size_t count_a = 0;
        std::size_t count_b = 0;
        for (std::size_t i = row; i < order; ++i)
        {
          count_a += (pattern[col][i] == 0 ? 0 : 1);
          count_b += (pattern[i][col] == 0 ? 0 : 1);
        }
        std::size_t count = (count_a - 1) * (count_b - 1);
        if (count < beta)
        {
          beta = count;
          max_row = col;
        }
      }
      // Swap row and max_row
      if (max_row != row)
      {
        for (std::size_t i = row; i < order; ++i)
          std::swap(pattern[row][i], pattern[max_row][i]);
        for (std::size_t i = row; i < order; ++i)
          std::swap(pattern[i][row], pattern[i][max_row]);
        std::swap(perm[row], perm[max_row]);
      }
      for (std::size_t col = row + 1; col < order; ++col)
        if (pattern[row][col])
          for (std::size_t i = row + 1; i < order; ++i)
            pattern[i][col] = pattern[i][row] || pattern[i][col];
    }
    return perm;
  }

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolver(
      const SparseMatrixPolicy& matrix,
      typename SparseMatrixPolicy::value_type initial_value)
      : LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>(
            matrix,
            initial_value,
            [](const SparseMatrixPolicy& m) -> LuDecompositionPolicy
            { return LuDecomposition::Create<SparseMatrixPolicy>(m); })
  {
  }

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::LinearSolver(
      const SparseMatrixPolicy& matrix,
      typename SparseMatrixPolicy::value_type initial_value,
      const std::function<LuDecompositionPolicy(const SparseMatrixPolicy&)> create_lu_decomp)
      : nLij_Lii_(),
        Lij_yj_(),
        nUij_Uii_(),
        Uij_xj_(),
        lu_decomp_(create_lu_decomp(matrix))
  {
    MICM_PROFILE_FUNCTION();

    auto lu = lu_decomp_.template GetLUMatrices<SparseMatrixPolicy>(matrix, initial_value);
    auto lower_matrix = std::move(lu.first);
    auto upper_matrix = std::move(lu.second);
    for (std::size_t i = 0; i < lower_matrix.NumRows(); ++i)
    {
      std::size_t nLij = 0;
      for (std::size_t j_id = lower_matrix.RowStartVector()[i]; j_id < lower_matrix.RowStartVector()[i + 1]; ++j_id)
      {
        std::size_t j = lower_matrix.RowIdsVector()[j_id];
        if (j >= i)
          break;
        Lij_yj_.push_back(std::make_pair(lower_matrix.VectorIndex(0, i, j), j));
        ++nLij;
      }
      // There must always be a non-zero element on the diagonal
      nLij_Lii_.push_back(std::make_pair(nLij, lower_matrix.VectorIndex(0, i, i)));
    }
    for (std::size_t i = upper_matrix.NumRows() - 1; i != static_cast<std::size_t>(-1); --i)
    {
      std::size_t nUij = 0;
      for (std::size_t j_id = upper_matrix.RowStartVector()[i]; j_id < upper_matrix.RowStartVector()[i + 1]; ++j_id)
      {
        std::size_t j = upper_matrix.RowIdsVector()[j_id];
        if (j <= i)
          continue;
        Uij_xj_.push_back(std::make_pair(upper_matrix.VectorIndex(0, i, j), j));
        ++nUij;
      }
      // There must always be a non-zero element on the diagonal
      nUij_Uii_.push_back(std::make_pair(nUij, upper_matrix.VectorIndex(0, i, i)));
    }
  };

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::Factor(
      const SparseMatrixPolicy& matrix,
      SparseMatrixPolicy& lower_matrix,
      SparseMatrixPolicy& upper_matrix)
  {
    MICM_PROFILE_FUNCTION();

    lu_decomp_.template Decompose<SparseMatrixPolicy>(matrix, lower_matrix, upper_matrix);
  }

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  inline void LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::Factor(
      const SparseMatrixPolicy& matrix,
      SparseMatrixPolicy& lower_matrix,
      SparseMatrixPolicy& upper_matrix,
      bool& is_singular)
  {
    MICM_PROFILE_FUNCTION();

    lu_decomp_.template Decompose<SparseMatrixPolicy>(matrix, lower_matrix, upper_matrix, is_singular);
  }

  template<class SparseMatrixPolicy, class LuDecompositionPolicy>
  template<class MatrixPolicy>
  requires(
      !VectorizableDense<MatrixPolicy> ||
      !VectorizableSparse<SparseMatrixPolicy>) inline void LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::
      Solve(
          MatrixPolicy& x,
          const SparseMatrixPolicy& lower_matrix,
          const SparseMatrixPolicy& upper_matrix) const
  {
    MICM_PROFILE_FUNCTION();

    for (std::size_t i_cell = 0; i_cell < x.NumRows(); ++i_cell)
    {
      auto x_cell = x[i_cell];
      const std::size_t lower_grid_offset = i_cell * lower_matrix.FlatBlockSize();
      const std::size_t upper_grid_offset = i_cell * upper_matrix.FlatBlockSize();
      auto& y_cell = x_cell;  // Alias x for consistency with equations, but to reuse memory
      
      // Forward Substitution
      {
        auto y_elem = y_cell.begin();
        auto Lij_yj = Lij_yj_.begin();
        for (auto& nLij_Lii : nLij_Lii_)
        {
          for (std::size_t i = 0; i < nLij_Lii.first; ++i)
          {
            *y_elem -= lower_matrix.AsVector()[lower_grid_offset + (*Lij_yj).first] * y_cell[(*Lij_yj).second];
            ++Lij_yj;
          }
          *(y_elem++) /= lower_matrix.AsVector()[lower_grid_offset + nLij_Lii.second];
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
            *x_elem -= upper_matrix.AsVector()[upper_grid_offset + (*Uij_xj).first] * x_cell[(*Uij_xj).second];
            ++Uij_xj;
          }

          *(x_elem) /= upper_matrix.AsVector()[upper_grid_offset + nUij_Uii.second];
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
  requires(VectorizableDense<MatrixPolicy>&&
               VectorizableSparse<SparseMatrixPolicy>) inline void LinearSolver<SparseMatrixPolicy, LuDecompositionPolicy>::
      Solve(
          MatrixPolicy& x,
          const SparseMatrixPolicy& lower_matrix,
          const SparseMatrixPolicy& upper_matrix) const
  {
    MICM_PROFILE_FUNCTION();
    const MatrixPolicy& b = x;
    const std::size_t n_cells = b.GroupVectorSize();
    // Loop over groups of blocks
    for (std::size_t i_group = 0; i_group < b.NumberOfGroups(); ++i_group)
    {
      auto b_group = std::next(b.AsVector().begin(), i_group * b.GroupSize());
      auto x_group = std::next(x.AsVector().begin(), i_group * x.GroupSize());
      auto L_group =
          std::next(lower_matrix.AsVector().begin(), i_group * lower_matrix.GroupSize(lower_matrix.FlatBlockSize()));
      auto U_group =
          std::next(upper_matrix.AsVector().begin(), i_group * upper_matrix.GroupSize(upper_matrix.FlatBlockSize()));
      {
        auto y_elem = x_group;
        auto Lij_yj = Lij_yj_.begin();
        for (auto& nLij_Lii : nLij_Lii_)
        {
          std::copy(b_group, b_group + n_cells, y_elem);
          b_group += n_cells;
          for (std::size_t i = 0; i < nLij_Lii.first; ++i)
          {
            std::size_t Lij_yj_first = (*Lij_yj).first;
            std::size_t Lij_yj_second_times_n_cells = (*Lij_yj).second * n_cells;
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              y_elem[i_cell] -= L_group[Lij_yj_first + i_cell] * x_group[Lij_yj_second_times_n_cells + i_cell];
            ++Lij_yj;
          }
          std::size_t nLij_Lii_second = nLij_Lii.second;
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            y_elem[i_cell] /= L_group[nLij_Lii_second + i_cell];
          y_elem += n_cells;
        }
      }
      {
        auto x_elem = std::next(x_group, x.GroupSize() - n_cells);
        auto Uij_xj = Uij_xj_.begin();
        for (auto& nUij_Uii : nUij_Uii_)
        {
          // x_elem starts out as y_elem from the previous loop
          for (std::size_t i = 0; i < nUij_Uii.first; ++i)
          {
            std::size_t Uij_xj_first = (*Uij_xj).first;
            std::size_t Uij_xj_second_times_n_cells = (*Uij_xj).second * n_cells;
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              x_elem[i_cell] -= U_group[Uij_xj_first + i_cell] * x_group[Uij_xj_second_times_n_cells + i_cell];
            ++Uij_xj;
          }
          std::size_t nUij_Uii_second = nUij_Uii.second;
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            x_elem[i_cell] /= U_group[nUij_Uii_second + i_cell];

          // don't iterate before the beginning of the vector
          std::size_t x_elem_distance = std::distance(x.AsVector().begin(), x_elem);
          x_elem -= std::min(n_cells, x_elem_distance);
        }
      }
    }
  }
}  // namespace micm