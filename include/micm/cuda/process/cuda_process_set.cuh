// Copyright (C) 2023-2025 National Center for Atmospheric Research
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

    /// This is the function that will copy the "jacobian_flat_id"
    ///   of class "ProcessSet" to the device, after the matrix
    ///   structure is known
    void CopyJacobiFlatId(ProcessSetParam& hoststruct, ProcessSetParam& devstruct);

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
  }  // namespace cuda
}  // namespace micm
