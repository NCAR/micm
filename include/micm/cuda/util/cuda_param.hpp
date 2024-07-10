// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <utility>

// To make the NormalizedError function works properly on GPU,
// make sure to choose the BLOCK_SIZE from [32, 64, 128, 256, 512, 1024]
const std::size_t BLOCK_SIZE = 32;

/// This struct holds the (1) pointer to, and (2) size of
///   each constatnt data member from the class "ProcessSet";
/// This struct could be allocated on the host or device;
struct ProcessSetParam
{
  std::size_t* number_of_reactants_ = nullptr;
  std::size_t* reactant_ids_ = nullptr;
  std::size_t* number_of_products_ = nullptr;
  std::size_t* product_ids_ = nullptr;
  double* yields_ = nullptr;
  std::size_t* jacobian_flat_ids_ = nullptr;
  std::size_t number_of_reactants_size_;
  std::size_t reactant_ids_size_;
  std::size_t number_of_products_size_;
  std::size_t product_ids_size_;
  std::size_t yields_size_;
  std::size_t jacobian_flat_ids_size_;
};

/// This struct holds the (1) pointer to, and (2) size of
///   each constatnt data member from the class "LuDecomposition";
/// This struct could be allocated on the host or device;
struct LuDecomposeParam
{
  bool* is_singular = nullptr;
  std::pair<std::size_t, std::size_t>* niLU_ = nullptr;
  char* do_aik_ = nullptr;
  std::size_t* aik_ = nullptr;
  std::pair<std::size_t, std::size_t>* uik_nkj_ = nullptr;
  std::pair<std::size_t, std::size_t>* lij_ujk_ = nullptr;
  char* do_aki_ = nullptr;
  std::size_t* aki_ = nullptr;
  std::pair<std::size_t, std::size_t>* lki_nkj_ = nullptr;
  std::pair<std::size_t, std::size_t>* lkj_uji_ = nullptr;
  std::size_t* uii_ = nullptr;
  std::size_t niLU_size_;
  std::size_t do_aik_size_;
  std::size_t aik_size_;
  std::size_t uik_nkj_size_;
  std::size_t lij_ujk_size_;
  std::size_t do_aki_size_;
  std::size_t aki_size_;
  std::size_t lki_nkj_size_;
  std::size_t lkj_uji_size_;
  std::size_t uii_size_;
};

/// This struct holds the (1) pointer to, and (2) size of
///   each constatnt data member from the class "LinearSolver";
/// This struct could be allocated on the host or device;
struct LinearSolverParam
{
  std::pair<std::size_t, std::size_t>* nLij_Lii_ = nullptr;
  std::pair<std::size_t, std::size_t>* Lij_yj_ = nullptr;
  std::pair<std::size_t, std::size_t>* nUij_Uii_ = nullptr;
  std::pair<std::size_t, std::size_t>* Uij_xj_ = nullptr;
  std::size_t nLij_Lii_size_;
  std::size_t Lij_yj_size_;
  std::size_t nUij_Uii_size_;
  std::size_t Uij_xj_size_;
};

/// This struct holds (1) pointer to, and (2) size of
///   data allocated on a device.
struct CudaMatrixParam
{
  double* d_data_ = nullptr;
  std::size_t number_of_elements_;
  std::size_t number_of_grid_cells_;
};

/// This struct holds (1) pointer to, and (2) size of
///   each constatnt data member from the class "CudaRosenbrockSolver";
/// This struct could be allocated on the host or device;
struct CudaRosenbrockSolverParam
{
  // for NormalizedError function
  double* errors_input_ = nullptr;
  double* errors_output_ = nullptr;
  double* absolute_tolerance_ = nullptr;
  std::size_t absolute_tolerance_size_;
  std::size_t errors_size_;
  // for AlphaMinusJacobian function
  std::size_t* jacobian_diagonal_elements_ = nullptr;
  std::size_t jacobian_diagonal_elements_size_;
};
