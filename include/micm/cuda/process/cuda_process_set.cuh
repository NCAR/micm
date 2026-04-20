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

    /// Copy pre-computed offset arrays for the optimized 2D kernel to device
    void CopyJacobianOffsets(
        const std::size_t* h_reactant_offsets,
        const std::size_t* h_flat_id_offsets,
        const std::size_t* h_yield_offsets,
        std::size_t offsets_size,
        ProcessSetParam& devstruct);

    /// Host function that calls the optimized 2D CUDA kernel
    ///   to form the negative Jacobian matrix (-J)
    void SubtractJacobianTermsKernel2DDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct);
    /// Narrow-cast and copy compact Jacobian arrays to device
    void CopyJacobianParamsCompact(
        const ProcessInfoParam* h_process_info,
        const std::size_t* h_reactant_ids,
        const std::size_t* h_flat_ids,
        std::size_t process_info_count,
        std::size_t reactant_ids_count,
        std::size_t flat_ids_count,
        ProcessSetParam& devstruct);

    /// Host function that calls the bandwidth-optimized compact kernel
    void SubtractJacobianTermsKernelCompactDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct);

    /// Host function for persistent kernel with parameterized block count.
    /// num_blocks controls cache-vs-parallelism tradeoff (smaller = more cache reuse).
    void SubtractJacobianTermsKernelPersistentDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct,
        int num_blocks);
  }  // namespace cuda
}  // namespace micm
