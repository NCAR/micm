// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <vector>

const size_t BLOCK_SIZE = 32;
// member data of class CudaProcessSet grouped in struct passing to kernel driver function
struct CudaProcessSetParam
{
  const size_t* number_of_reactants_;
  const size_t* reactant_ids_;
  size_t reactant_ids_size_;
  const size_t* number_of_products_;
  const size_t* product_ids_;
  size_t product_ids_size_;
  const double* yields_;
  size_t yields_size_;
  const size_t* jacobian_flat_ids_;
  size_t jacobian_flat_ids_size_;
};

struct CudaSolverParam
{
  const std::pair<size_t, size_t>* niLU_;
  size_t niLU_size_;
  const char* do_aik_;
  size_t do_aik_size_;
  const size_t* aik_;
  size_t aik_size_;
  const std::pair<size_t, size_t>* uik_nkj_;
  size_t uik_nkj_size_;
  const std::pair<size_t, size_t>* lij_ujk_;
  size_t lij_ujk_size_;
  const char* do_aki_;
  size_t do_aki_size_;
  const size_t* aki_;
  size_t aki_size_;
  const std::pair<size_t, size_t>* lki_nkj_;
  size_t lki_nkj_size_;
  const std::pair<size_t, size_t>* lkj_uji_;
  size_t lkj_uji_size_;
  const size_t* uii_;
  size_t uii_size_;
};

struct CudaLinearSolverParam
{
  const std::pair<size_t, size_t>* nLij_Lii_;
  size_t nLij_Lii_size_;
  const std::pair<size_t, size_t>* Lij_yj_;
  size_t Lij_yj_size_;
  const std::pair<size_t, size_t>* nUij_Uii_;
  size_t nUij_Uii_size_;
  const std::pair<size_t, size_t>* Uij_xj_;
  size_t Uij_xj_size_;
};

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

/// The device struct that holds the constatnt data members 
/// from the "CudaLuDecomposition" class type;
struct LuDecomposeConstDevice
{
  std::pair<size_t, size_t>* d_niLU_;
  char* d_do_aik_;
  size_t* d_aik_;
  std::pair<size_t, size_t>* d_uik_nkj_;
  std::pair<size_t, size_t>* d_lij_ujk_;
  char* d_do_aki_;
  size_t* d_aki_;
  std::pair<size_t, size_t>* d_lki_nkj_;
  std::pair<size_t, size_t>* d_lkj_uji_;
  size_t* d_uii_;
};
