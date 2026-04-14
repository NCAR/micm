// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

/// @file host_device.hpp
/// @brief Portability macro for functions that must compile on both CPU and GPU.
///
/// MICM_HOST_DEVICE expands to __host__ __device__ when compiling with nvcc,
/// and to nothing otherwise. Apply it to any free function whose body is
/// shared between CPU and CUDA kernels (e.g. CalculateBatch overloads in
/// rate_constant_functions.hpp).

#ifdef __CUDACC__
#  define MICM_HOST_DEVICE __host__ __device__
#else
#  define MICM_HOST_DEVICE
#endif
