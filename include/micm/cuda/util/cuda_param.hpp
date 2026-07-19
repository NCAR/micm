// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/util/types.hpp>

#include <cstddef>
#include <cstdint>
#include <utility>

// To make the NormalizedError function works properly on GPU,
// make sure to choose the BLOCK_SIZE from [32, 64, 128, 256, 512, 1024]
const micm::Index BLOCK_SIZE = 512;

/// This struct holds information about a process for the Jacobian calculation
struct ProcessInfoParam
{
  micm::Index process_id_;
  micm::Index independent_id_;
  micm::Index number_of_dependent_reactants_;
  micm::Index number_of_products_;
};

/// This struct holds the (1) pointer to, and (2) size of
///   each constatnt data member from the class "ProcessSet";
/// This struct could be allocated on the host or device;
struct ProcessSetParam
{
  micm::Index* number_of_reactants_ = nullptr;
  micm::Index* reactant_ids_ = nullptr;
  micm::Index* number_of_products_ = nullptr;
  micm::Index* product_ids_ = nullptr;
  micm::Real* yields_ = nullptr;
  ProcessInfoParam* jacobian_process_info_ = nullptr;
  micm::Index* jacobian_reactant_ids_ = nullptr;
  micm::Index* jacobian_product_ids_ = nullptr;
  micm::Real* jacobian_yields_ = nullptr;
  micm::Index* jacobian_flat_ids_ = nullptr;
  uint8_t* is_algebraic_variable_ = nullptr;
  micm::Index number_of_reactants_size_;
  micm::Index reactant_ids_size_;
  micm::Index number_of_products_size_;
  micm::Index product_ids_size_;
  micm::Index yields_size_;
  micm::Index jacobian_process_info_size_;
  micm::Index jacobian_reactant_ids_size_;
  micm::Index jacobian_product_ids_size_;
  micm::Index jacobian_yields_size_;
  micm::Index jacobian_flat_ids_size_;
  micm::Index algebraic_variable_size_;
};

/// This struct holds the (1) pointer to, and (2) size of
///   each constant data member from the class "LuDecompositionMozartInPlace";
/// This struct could be allocated on the host or device;
struct LuDecomposeMozartInPlaceParam
{
  std::tuple<micm::Index, micm::Index, micm::Index>* aii_nji_nki_ = nullptr;
  micm::Index* aji_ = nullptr;
  std::pair<micm::Index, micm::Index>* aik_njk_ = nullptr;
  std::pair<micm::Index, micm::Index>* ajk_aji_ = nullptr;
  micm::Index aii_nji_nki_size_;
  micm::Index aji_size_;
  micm::Index aik_njk_size_;
  micm::Index ajk_aji_size_;
  micm::Index number_of_non_zeros_;
};

/// Alias for the default LU decomposition parameter struct
using LuDecomposeParam = LuDecomposeMozartInPlaceParam;

/// This struct holds the (1) pointer to, and (2) size of
///   each constatnt data member from the class "LinearSolverInPlace";
/// This struct could be allocated on the host or device;
struct LinearSolverInPlaceParam
{
  micm::Index* nLij_ = nullptr;
  std::pair<micm::Index, micm::Index>* Lij_yj_ = nullptr;
  std::pair<micm::Index, micm::Index>* nUij_Uii_ = nullptr;
  std::pair<micm::Index, micm::Index>* Uij_xj_ = nullptr;
  micm::Index nLij_size_;
  micm::Index Lij_yj_size_;
  micm::Index nUij_Uii_size_;
  micm::Index Uij_xj_size_;
  micm::Index number_of_non_zeros_;
};

/// This struct holds (1) pointer to, and (2) size of
///   data allocated on a device.
struct CudaMatrixParam
{
  micm::Real* d_data_ = nullptr;
  micm::Index number_of_elements_;
  micm::Index number_of_grid_cells_;
  micm::Index vector_length_;
};

struct CudaErrorParam
{
  micm::Real* errors_input_ = nullptr;
  micm::Real* errors_output_ = nullptr;
  micm::Index errors_size_;
};

struct CudaJacobianDiagonalElementsParam
{
  micm::Index* data_ = nullptr;
  micm::Index size_;
};