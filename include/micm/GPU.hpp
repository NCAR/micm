// Copyright (C) 2024-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/cuda/process/cuda_process_set.cuh>
#include <micm/cuda/process/cuda_process_set.hpp>
#include <micm/cuda/solver/cuda_linear_solver_in_place.cuh>
#include <micm/cuda/solver/cuda_linear_solver_in_place.hpp>
#include <micm/cuda/solver/cuda_lu_decomposition.hpp>
#include <micm/cuda/solver/cuda_lu_decomposition_mozart_in_place.cuh>
#include <micm/cuda/solver/cuda_lu_decomposition_mozart_in_place.hpp>
#include <micm/cuda/solver/cuda_rosenbrock.cuh>
#include <micm/cuda/solver/cuda_rosenbrock.hpp>
#include <micm/cuda/solver/cuda_solver_builder.hpp>
#include <micm/cuda/solver/cuda_solver_parameters.hpp>
#include <micm/cuda/solver/cuda_state.hpp>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_matrix.cuh>
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/cuda/util/cuda_util.cuh>

namespace micm
{
  template<std::size_t L>
  using GpuBuilder = micm::CudaSolverBuilderInPlace<micm::CudaRosenbrockSolverParameters, L>;
}