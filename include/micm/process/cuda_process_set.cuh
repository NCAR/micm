// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <micm/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    /// This is the function that will copy the constant data
    ///   members of class "ProcessSet" to the device, except
    ///   for the "jacobian_flat_id" because it is unknown now
    ProcessSetParam CopyConstData(ProcessSetParam& hoststruct);

    /// This is the function that will copy the "jacobian_flat_id"
    ///   of class "ProcessSet" to the device, after the matrix
    ///   structure is known
    void CopyJacobiFlatId(ProcessSetParam& hoststruct, ProcessSetParam& devstruct);

    /// This is the host function that will call the CUDA kernel
    ///   to calculate the forcing terms
    void AddForcingTermsKernelDriver(
        const CudaVectorMatrixParam& rate_constants_param,
        const CudaVectorMatrixParam& state_variables_param,
        CudaVectorMatrixParam& forcing_param,
        const ProcessSetParam& devstruct);

    /// This is the host function that will call the CUDA kernel
    ///   to form the negative Jacobian matrix (-J)
    std::chrono::nanoseconds SubtractJacobianTermsKernelDriver(
        CudaMatrixParam& matrixParam, 
        CudaSparseMatrixParam& sparseMatrix, 
        const ProcessSetParam& devstruct);
  }  // namespace cuda
}  // namespace micm
