// Copyright (C) 2023 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

namespace micm
{

  /// @brief A general-use sparse-matrix linear solver
  template<typename T, template<class> class SparseMatrixPolicy>
  class LinearSolver
  {
    // Parameters needed to calculate L (U x) = b
    //
    // The calculation is split into calculation of L y = b where y = U x:
    //
    // y_1 = b_1 / L_11
    // y_i = 1 / L_ii * [ b_i - sum( j = 1...i-1 ){ L_ij * y_j } ]  i = 2...N
    //
    // ... and then U x = y:
    //
    // x_N = y_N / U_NN
    // x_i = 1 / U_ii * [ y_i - sum( j = i+1...N ){ U_ij * x_j } ] i = N-1...1

    // Number of non-zero elements (excluding the diagonal) and the index of the diagonal
    // element for each row in L
    std::vector<std::pair<std::size_t, std::size_t>> nLij_Lii_;
    // Indices of non-zero combinations of L_ij and y_j
    std::vector<std::pair<std::size_t, std::size_t>> Lij_yj_;
    // Number of non-zero elements (exluding the diagonal) and the index of the diagonal
    // element for each row in U (in reverse order)
    std::vector<std::pair<std::size_t, std::size_t>> nUij_Uii_;
    // Indices of non-zero combinations of U_ij and x_j
    std::vector<std::pair<std::size_t, std::size_t>> Uij_xj_;

    LuDecomposition lu_decomp_;
    SparseMatrixPolicy<T> lower_matrix_;
    SparseMatrixPolicy<T> upper_matrix_;

   public:
    /// @brief default constructor
    LinearSolver() = default;

    /// @brief Constructs a linear solver for the sparsity structure of the given matrix
    /// @param matrix Sparse matrix
    LinearSolver(const SparseMatrixPolicy<T>& matrix, T initial_value);

    /// @brief Decompose the matrix into upper and lower triangular matrices
    void Factor(SparseMatrixPolicy<T>& matrix);

    /// @brief Solve for x in Ax = b
    template<template<class> class MatrixPolicy>
      requires(!VectorizableDense<MatrixPolicy<T>> || !VectorizableSparse<SparseMatrixPolicy<T>>)
    void Solve(const MatrixPolicy<T>& b, MatrixPolicy<T>& x);
    template<template<class> class MatrixPolicy>
      requires(VectorizableDense<MatrixPolicy<T>> && VectorizableSparse<SparseMatrixPolicy<T>>)
    void Solve(const MatrixPolicy<T>& b, MatrixPolicy<T>& x);
  };

  template<typename T, template<class> class SparseMatrixPolicy>
  inline LinearSolver<T, SparseMatrixPolicy>::LinearSolver(const SparseMatrixPolicy<T>& matrix, T initial_value)
      : nLij_Lii_(),
        Lij_yj_(),
        nUij_Uii_(),
        Uij_xj_(),
        lu_decomp_(matrix)
  {
    auto lu = lu_decomp_.GetLUMatrices(matrix, initial_value);
    lower_matrix_ = std::move(lu.first);
    upper_matrix_ = std::move(lu.second);
    for (std::size_t i = 0; i < lower_matrix_[0].size(); ++i)
    {
      std::size_t nLij = 0;
      for (std::size_t j_id = lower_matrix_.RowStartVector()[i]; j_id < lower_matrix_.RowStartVector()[i + 1]; ++j_id)
      {
        std::size_t j = lower_matrix_.RowIdsVector()[j_id];
        if (j >= i)
          break;
        Lij_yj_.push_back(std::make_pair(lower_matrix_.VectorIndex(0, i, j), j));
        ++nLij;
      }
      // There must always be a non-zero element on the diagonal
      nLij_Lii_.push_back(std::make_pair(nLij, lower_matrix_.VectorIndex(0, i, i)));
    }
    for (std::size_t i = upper_matrix_[0].size() - 1; i != static_cast<std::size_t>(-1); --i)
    {
      std::size_t nUij = 0;
      for (std::size_t j_id = upper_matrix_.RowStartVector()[i]; j_id < upper_matrix_.RowStartVector()[i + 1]; ++j_id)
      {
        std::size_t j = upper_matrix_.RowIdsVector()[j_id];
        if (j <= i)
          continue;
        Uij_xj_.push_back(std::make_pair(upper_matrix_.VectorIndex(0, i, j), j));
        ++nUij;
      }
      // There must always be a non-zero element on the diagonal
      nUij_Uii_.push_back(std::make_pair(nUij, upper_matrix_.VectorIndex(0, i, i)));
    }
  };

  template<typename T, template<class> class SparseMatrixPolicy>
  inline void LinearSolver<T, SparseMatrixPolicy>::Factor(SparseMatrixPolicy<T>& matrix)
  {
    lu_decomp_.Decompose<T, SparseMatrixPolicy>(matrix, lower_matrix_, upper_matrix_);
  }

  template<typename T, template<class> class SparseMatrixPolicy>
  template<template<class> class MatrixPolicy>
    requires(!VectorizableDense<MatrixPolicy<T>> || !VectorizableSparse<SparseMatrixPolicy<T>>)
  inline void LinearSolver<T, SparseMatrixPolicy>::Solve(const MatrixPolicy<T>& b, MatrixPolicy<T>& x)
  {
    for (std::size_t i_cell = 0; i_cell < b.size(); ++i_cell)
    {
      auto b_cell = b[i_cell];
      auto x_cell = x[i_cell];
      const std::size_t lower_grid_offset = i_cell * lower_matrix_.FlatBlockSize();
      const std::size_t upper_grid_offset = i_cell * upper_matrix_.FlatBlockSize();
      auto& y_cell = x_cell;  // Alias x for consistency with equations, but to reuse memory
      {
        auto b_elem = b_cell.begin();
        auto y_elem = y_cell.begin();
        auto Lij_yj = Lij_yj_.begin();
        for (auto& nLij_Lii : nLij_Lii_)
        {
          *y_elem = *(b_elem++);
          for (std::size_t i = 0; i < nLij_Lii.first; ++i)
          {
            *y_elem -= lower_matrix_.AsVector()[lower_grid_offset + (*Lij_yj).first] * y_cell[(*Lij_yj).second];
            ++Lij_yj;
          }
          *(y_elem++) /= lower_matrix_.AsVector()[lower_grid_offset + nLij_Lii.second];
        }
      }
      {
        auto y_elem = std::next(y_cell.end(), -1);
        auto x_elem = std::next(x_cell.end(), -1);
        auto Uij_xj = Uij_xj_.begin();
        for (auto& nUij_Uii : nUij_Uii_)
        {
          *x_elem = *(y_elem--);
          for (std::size_t i = 0; i < nUij_Uii.first; ++i)
          {
            *x_elem -= upper_matrix_.AsVector()[upper_grid_offset + (*Uij_xj).first] * x_cell[(*Uij_xj).second];
            ++Uij_xj;
          }
          *(x_elem--) /= upper_matrix_.AsVector()[upper_grid_offset + nUij_Uii.second];
        }
      }
    }
  }

  template<typename T, template<class> class SparseMatrixPolicy>
  template<template<class> class MatrixPolicy>
    requires(VectorizableDense<MatrixPolicy<T>> && VectorizableSparse<SparseMatrixPolicy<T>>)
  inline void LinearSolver<T, SparseMatrixPolicy>::Solve(const MatrixPolicy<T>& b, MatrixPolicy<T>& x)
  {
    const std::size_t n_cells = b.GroupVectorSize();
    // Loop over groups of blocks
    for(std::size_t i_group = 0; i_group < b.NumberOfGroups(); ++i_group)
    {
      auto b_group = std::next(b.AsVector().begin(), i_group * b.GroupSize());
      auto x_group = std::next(x.AsVector().begin(), i_group * x.GroupSize());
      auto L_group = std::next(lower_matrix_.AsVector().begin(), i_group * lower_matrix_.GroupSize(lower_matrix_.FlatBlockSize()));
      auto U_group = std::next(upper_matrix_.AsVector().begin(), i_group * upper_matrix_.GroupSize(upper_matrix_.FlatBlockSize()));
      auto y_group = x_group; // Alias x for consistency with equations, but to reuse memory
      {
        auto b_elem = b_group;
        auto y_elem = y_group;
        auto Lij_yj = Lij_yj_.begin();
        for (auto& nLij_Lii : nLij_Lii_)
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            y_elem[i_cell] = b_elem[i_cell];
          b_elem += n_cells;          
          for (std::size_t i = 0; i < nLij_Lii.first; ++i)
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              y_elem[i_cell] -= L_group[(*Lij_yj).first + i_cell] * y_group[(*Lij_yj).second * n_cells + i_cell];
            ++Lij_yj;
          }
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            y_elem[i_cell] /= L_group[nLij_Lii.second + i_cell];
          y_elem += n_cells;
        }
      }
      {
        auto y_elem = std::next(y_group, x.GroupSize() - n_cells);
        auto x_elem = std::next(x_group, x.GroupSize() - n_cells);
        auto Uij_xj = Uij_xj_.begin();
        for (auto& nUij_Uii : nUij_Uii_)
        {
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            x_elem[i_cell] = y_elem[i_cell];
          y_elem -= n_cells;
          for (std::size_t i = 0; i < nUij_Uii.first; ++i)
          {
            for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
              x_elem[i_cell] -= U_group[(*Uij_xj).first + i_cell] * x_group[(*Uij_xj).second * n_cells + i_cell];
            ++Uij_xj;
          }
          for (std::size_t i_cell = 0; i_cell < n_cells; ++i_cell)
            x_elem[i_cell] /= U_group[nUij_Uii.second + i_cell];
          x_elem -= n_cells;
        }
      }
    }
  }

}  // namespace micm