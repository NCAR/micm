// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_param.hpp>

namespace micm
{
  namespace cuda
  {
    /// @brief Copy rate constant descriptors to device memory
    /// @param hoststruct Host-side struct with rate constant data pointers
    /// @return Device-side struct with allocated and copied device memory
    ProcessParam CopyProcessConstData(ProcessParam& hoststruct);

    /// @brief Free device memory for rate constant descriptors
    /// @param devstruct Device-side struct to free
    void FreeProcessConstData(ProcessParam& devstruct);

    /// @brief Launch CUDA kernel to calculate rate constants
    /// @param d_temperature Device array of temperatures [K] (one per grid cell)
    /// @param d_pressure Device array of pressures [Pa] (one per grid cell)
    /// @param d_air_density Device array of air densities [mol m-3] (one per grid cell)
    /// @param custom_rate_params Custom rate parameters matrix on device
    /// @param d_fixed_reactants Device array of pre-computed parameterized reactant factors
    /// @param rate_constants_param Output rate constants matrix on device
    /// @param devstruct Device-side struct with rate constant descriptors
    void CalculateRateConstantsKernelDriver(
        const double* d_temperature,
        const double* d_pressure,
        const double* d_air_density,
        const CudaMatrixParam& custom_rate_params,
        const double* d_fixed_reactants,
        CudaMatrixParam& rate_constants_param,
        const ProcessParam& devstruct);
  }  // namespace cuda
}  // namespace micm
