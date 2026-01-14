// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <utility>

// To make the NormalizedError function works properly on GPU,
// make sure to choose the BLOCK_SIZE from [32, 64, 128, 256, 512, 1024]
const std::size_t BLOCK_SIZE = 512;

/// This struct holds information about a process for the Jacobian calculation
struct ProcessInfoParam
{
  std::size_t process_id_;
  std::size_t independent_id_;
  std::size_t number_of_dependent_reactants_;
  std::size_t number_of_products_;
};

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
  ProcessInfoParam* jacobian_process_info_ = nullptr;
  std::size_t* jacobian_reactant_ids_ = nullptr;
  std::size_t* jacobian_product_ids_ = nullptr;
  double* jacobian_yields_ = nullptr;
  std::size_t* jacobian_flat_ids_ = nullptr;
  std::size_t number_of_reactants_size_;
  std::size_t reactant_ids_size_;
  std::size_t number_of_products_size_;
  std::size_t product_ids_size_;
  std::size_t yields_size_;
  std::size_t jacobian_process_info_size_;
  std::size_t jacobian_reactant_ids_size_;
  std::size_t jacobian_product_ids_size_;
  std::size_t jacobian_yields_size_;
  std::size_t jacobian_flat_ids_size_;
};

/// This struct holds the (1) pointer to, and (2) size of
///   each constant data member from the class "LuDecompositionMozartInPlace";
/// This struct could be allocated on the host or device;
struct LuDecomposeMozartInPlaceParam
{
  std::tuple<std::size_t, std::size_t, std::size_t>* aii_nji_nki_ = nullptr;
  std::size_t* aji_ = nullptr;
  std::pair<std::size_t, std::size_t>* aik_njk_ = nullptr;
  std::pair<std::size_t, std::size_t>* ajk_aji_ = nullptr;
  std::size_t aii_nji_nki_size_;
  std::size_t aji_size_;
  std::size_t aik_njk_size_;
  std::size_t ajk_aji_size_;
  std::size_t number_of_non_zeros_;
};

/// Alias for the default LU decomposition parameter struct
using LuDecomposeParam = LuDecomposeMozartInPlaceParam;

/// This struct holds the (1) pointer to, and (2) size of
///   each constatnt data member from the class "LinearSolverInPlace";
/// This struct could be allocated on the host or device;
struct LinearSolverInPlaceParam
{
  std::size_t* nLij_ = nullptr;
  std::pair<std::size_t, std::size_t>* Lij_yj_ = nullptr;
  std::pair<std::size_t, std::size_t>* nUij_Uii_ = nullptr;
  std::pair<std::size_t, std::size_t>* Uij_xj_ = nullptr;
  std::size_t nLij_size_;
  std::size_t Lij_yj_size_;
  std::size_t nUij_Uii_size_;
  std::size_t Uij_xj_size_;
  std::size_t number_of_non_zeros_;
};

/// This struct holds (1) pointer to, and (2) size of
///   data allocated on a device.
struct CudaMatrixParam
{
  double* d_data_ = nullptr;
  std::size_t number_of_elements_;
  std::size_t number_of_grid_cells_;
  std::size_t vector_length_;
};

struct CudaErrorParam
{
  double* errors_input_ = nullptr;
  double* errors_output_ = nullptr;
  std::size_t errors_size_;
};

struct CudaJacobianDiagonalElementsParam
{
  std::size_t* data_ = nullptr;
  std::size_t size_;
};