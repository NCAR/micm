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

    /// Host function that calls the shared-memory CUDA kernel to form the negative Jacobian matrix (-J).
    /// Layout: 1 block = 1 grid cell; threads in a block partition the i_proc loop and
    /// accumulate contributions into a per-cell shared-memory buffer using shared-memory
    /// atomicAdd (NO atomics on global d_jacobian). Requires the three jacobian_*_offsets_
    /// arrays in devstruct to be populated.
    void SubtractJacobianTermsKernelSharedDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct);

    /// Host function that calls the deterministic-reduction shared-memory CUDA kernel
    /// to form the negative Jacobian matrix (-J). Strides i_proc across threads, builds
    /// per-thread register-resident contribution lists, and merges them into a per-cell
    /// shared-memory accumulator via a chunked warp-shuffle tree reduction. NO atomics
    /// (neither global nor shared). Reduction order is fixed across runs.
    void SubtractJacobianTermsKernelSharedReduceDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct);
  }  // namespace cuda
}  // namespace micm
