// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    /// This is the function that will copy the constant data
    ///   members of class "ProcessSet" to the device, except
    ///   for the "jacobian_flat_id" because it is unknown now
    ProcessSetParam CopyConstData(ProcessSetParam& hoststruct);

    /// This is the function that will delete the constant data
    ///   members of class "CudaProcessSet" on the device
    void FreeConstData(ProcessSetParam& devstruct);

    /// This is the function that will copy the information needed
    /// to calculated Jacobian terms to the device
    void CopyJacobianParams(ProcessSetParam& hoststruct, ProcessSetParam& devstruct);

    /// This is the host function that will call the CUDA kernel
    ///   to calculate the forcing terms
    void AddForcingTermsKernelDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& forcing_param,
        const ProcessSetParam& devstruct);

    /// This is the host function that will call the CUDA kernel
    ///   to form the negative Jacobian matrix (-J)
    void SubtractJacobianTermsKernelDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct);

    /// Host function that calls the unrolled CUDA kernel to form the negative Jacobian matrix (-J).
    /// Same data layout as the original kernel, but with #pragma unroll hints on all loops
    /// (outer process_info loop unrolled by 4; inner loops fully unrolled when bounds are known).
    void SubtractJacobianTermsKernelUnrolledDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct);

    /// Host function that calls the 2D CUDA kernel to form the negative Jacobian matrix (-J).
    /// Launches a 2D grid: X dimension over grid cells, Y dimension over process_infos.
    /// Each thread handles one (cell_id, i_proc) pair independently. Contributions from
    /// different i_proc threads targeting the same Jacobian entry are merged via global
    /// atomicAdd. Requires the three jacobian_*_offsets_ arrays in devstruct to be populated.
    void SubtractJacobianTermsKernel2DDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct);

    /// Single-pass gather kernel driver. Parallelizes on (grid cell × unique Jacobian entry).
    /// Each thread owns one Jacobian entry, loops over its CSR contributors, recomputes
    /// d_rate_d_ind per contributor, accumulates, and writes once — no atomics.
    /// Requires jac_gather_* arrays in devstruct (populated by SetJacobianFlatIds).
    void SubtractJacobianTermsGatherKernelDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct);

    /// Two-pass gather driver.
    /// Pass 1 (1D grid, same shape as the original scatter kernel): compute d_rate_d_ind for
    ///   every process_info and store in a scratch buffer sized N_groups×N_proc×VL.
    /// Pass 2 (2D grid over cells × unique Jacobian entries): read from scratch via the CSR,
    ///   accumulate weighted contributions, and write each Jacobian entry once — no atomics.
    /// Allocates and frees the scratch buffer per call (cudaMallocAsync/cudaFreeAsync on
    /// stream 0). Requires jac_gather_* arrays in devstruct (populated by SetJacobianFlatIds).
    void SubtractJacobianTermsTwoPassGatherDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct);
  }  // namespace cuda
}  // namespace micm
