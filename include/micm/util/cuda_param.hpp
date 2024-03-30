// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

// To make the NormalizedError function works properly on GPU,
// make sure to choose the BLOCK_SIZE from [32, 64, 128, 256, 512, 1024]
const size_t BLOCK_SIZE = 32;

// different matrix data grouped in struct passing to kernel driver function
struct CudaMatrixParam
{
  const double* rate_constants_;
  const double* state_variables_;
  double* forcing_;
  const double* b_;
  double* x_;
  size_t x_size_;
  size_t b_size_;
  size_t n_grids_;
  size_t n_reactions_;
  size_t n_species_;
  size_t b_column_counts_;
  size_t x_column_counts_;
};

// sparseMatrix data grouped in struct passing to kernel driver function
struct CudaSparseMatrixParam
{
  double* jacobian_;
  size_t jacobian_size_;
  const double* A_;
  size_t A_size_;
  double* L_;
  size_t L_size_;
  double* U_;
  size_t U_size_;
  size_t n_grids_;
  const double* lower_matrix_;
  size_t lower_matrix_size_;
  const double* upper_matrix_;
  size_t upper_matrix_size_;
};

/// This struct holds the (1) pointer to, and (2) size of
///   each constatnt data member from the class "ProcessSet";
/// This struct could be allocated on the host or device;
struct ProcessSetParam
{
  size_t* number_of_reactants_;
  size_t* reactant_ids_;
  size_t* number_of_products_;
  size_t* product_ids_;
  double* yields_;
  size_t* jacobian_flat_ids_;
  size_t number_of_reactants_size_;
  size_t reactant_ids_size_;
  size_t number_of_products_size_;
  size_t product_ids_size_;
  size_t yields_size_;
  size_t jacobian_flat_ids_size_;
};

/// This struct holds the (1) pointer to, and (2) size of
///   each constatnt data member from the class "LuDecomposition";
/// This struct could be allocated on the host or device;
struct LuDecomposeParam
{
  std::pair<size_t, size_t>* niLU_;
  char* do_aik_;
  size_t* aik_;
  std::pair<size_t, size_t>* uik_nkj_;
  std::pair<size_t, size_t>* lij_ujk_;
  char* do_aki_;
  size_t* aki_;
  std::pair<size_t, size_t>* lki_nkj_;
  std::pair<size_t, size_t>* lkj_uji_;
  size_t* uii_;
  size_t niLU_size_;
  size_t do_aik_size_;
  size_t aik_size_;
  size_t uik_nkj_size_;
  size_t lij_ujk_size_;
  size_t do_aki_size_;
  size_t aki_size_;
  size_t lki_nkj_size_;
  size_t lkj_uji_size_;
  size_t uii_size_;
};

/// This struct holds the (1) pointer to, and (2) size of
///   each constatnt data member from the class "LinearSolver";
/// This struct could be allocated on the host or device;
struct LinearSolverParam
{
  std::pair<size_t, size_t>* nLij_Lii_;
  std::pair<size_t, size_t>* Lij_yj_;
  std::pair<size_t, size_t>* nUij_Uii_;
  std::pair<size_t, size_t>* Uij_xj_;
  size_t nLij_Lii_size_;
  size_t Lij_yj_size_;
  size_t nUij_Uii_size_;
  size_t Uij_xj_size_;
};

/// This struct holds (1) pointer to, and (2) size of
///   data allocated on a device.
struct CudaVectorMatrixParam
{
  double* d_data_;
  size_t number_of_elements_;
  size_t number_of_grid_cells_;
};

/// This struct holds (1) pointer to, and (2) size of
///   each constatnt data member from the class "CudaRosenbrockSolver";
/// This struct could be allocated on the host or device;
struct CudaRosenbrockSolverParam
{
  size_t num_grid_cells_;
  // for NormalizedError function
  double* errors_input_;
  double* errors_output_;
  size_t errors_size_;
  // for AlphaMinusJacobian function
  size_t* jacobian_diagonal_elements_;
  size_t jacobian_diagonal_elements_size_;
};
