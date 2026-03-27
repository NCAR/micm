// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>

#define _USE_MATH_DEFINES
#include <cmath>

namespace micm
{
  namespace cuda
  {
    // ========================================================================
    // Device functions for each rate constant type
    // ========================================================================

    /// @brief Arrhenius: k = A * exp(C/T) * (T/D)^B * (1 + E*P)
    /// params: [A, B, C, D, E]
    __device__ double CalculateArrhenius(double T, double P, const double* params)
    {
      return params[0] * exp(params[2] / T) * pow(T / params[3], params[1]) * (1.0 + params[4] * P);
    }

    /// @brief Troe falloff: k0*[M]/(1+k0*[M]/kinf) * Fc^(N/(N+log10(k0*[M]/kinf)^2))
    /// params: [k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N]
    __device__ double CalculateTroe(double T, double air_density, const double* params)
    {
      double k0 = params[0] * exp(params[2] / T) * pow(T / 300.0, params[1]);
      double kinf = params[3] * exp(params[5] / T) * pow(T / 300.0, params[4]);
      double ratio = k0 * air_density / kinf;
      return k0 * air_density / (1.0 + ratio) * pow(params[6], params[7] / (params[7] + pow(log10(ratio), 2)));
    }

    /// @brief Tunneling: k = A * exp(-B/T + C/T^3)
    /// params: [A, B, C]
    __device__ double CalculateTunneling(double T, const double* params)
    {
      return params[0] * exp(-params[1] / T + params[2] / (T * T * T));
    }

    /// @brief Branched rate constant
    /// params: [X, Y, k0, z, branch(0=Alkoxy,1=Nitrate)]
    /// k0 and z are precomputed on host during construction
    __device__ double CalculateBranched(double T, double air_density, const double* params)
    {
      double pre = params[0] * exp(-params[1] / T);
      double k0 = params[2];
      double z = params[3];
      int branch = static_cast<int>(params[4]);

      // Calculate A(T,[M])
      double a = k0 * air_density;
      double b = 0.43 * pow(T / 298.0, -8.0);
      double Atmn = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2)));

      if (branch == 0)  // Alkoxy
        return pre * (z / (z + Atmn));
      else  // Nitrate
        return pre * (Atmn / (Atmn + z));
    }

    /// @brief Ternary Chemical Activation: k0/(1+k0*[M]/kinf) * Fc^(N/(N+log10(k0*[M]/kinf)^2))
    /// params: [k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N]
    __device__ double CalculateTernaryChemicalActivation(double T, double air_density, const double* params)
    {
      double k0 = params[0] * exp(params[2] / T) * pow(T / 300.0, params[1]);
      double kinf = params[3] * exp(params[5] / T) * pow(T / 300.0, params[4]);
      double ratio = k0 * air_density / kinf;
      return k0 / (1.0 + ratio) * pow(params[6], params[7] / (params[7] + pow(log10(ratio), 2)));
    }

    /// @brief Reversible: k = A * exp(C/T) * k_r
    /// params: [A, C, k_r]
    __device__ double CalculateReversible(double T, const double* params)
    {
      double K_eq = params[0] * exp(params[1] / T);
      return K_eq * params[2];
    }

    /// @brief Surface: k = 4*N*pi*r^2 / (r/D + 4/(v_mean*gamma))
    /// params: [diffusion_coeff, mean_free_speed_factor, reaction_probability]
    /// custom_params: [radius, number_concentration]
    __device__ double CalculateSurface(double T, const double* params, const double* custom_params)
    {
      double mean_free_speed = sqrt(params[1] * T);
      double radius = custom_params[0];
      double number = custom_params[1];
      return 4.0 * number * M_PI * radius * radius /
             (radius / params[0] + 4.0 / (mean_free_speed * params[2]));
    }

    /// @brief UserDefined: k = custom_params[0] * scaling_factor
    /// params: [scaling_factor]
    __device__ double CalculateUserDefined(const double* params, const double* custom_params)
    {
      return custom_params[0] * params[0];
    }

    /// @brief Dispatch to the appropriate rate constant calculation
    __device__ double CalculateRateConstant(
        const CudaRateConstantData& rc,
        double T,
        double P,
        double air_density,
        const double* custom_params)
    {
      switch (rc.type_)
      {
        case CudaRateConstantType::Arrhenius: return CalculateArrhenius(T, P, rc.params_);
        case CudaRateConstantType::Troe: return CalculateTroe(T, air_density, rc.params_);
        case CudaRateConstantType::Tunneling: return CalculateTunneling(T, rc.params_);
        case CudaRateConstantType::Branched: return CalculateBranched(T, air_density, rc.params_);
        case CudaRateConstantType::TernaryChemicalActivation:
          return CalculateTernaryChemicalActivation(T, air_density, rc.params_);
        case CudaRateConstantType::Reversible: return CalculateReversible(T, rc.params_);
        case CudaRateConstantType::Surface: return CalculateSurface(T, rc.params_, custom_params);
        case CudaRateConstantType::UserDefined: return CalculateUserDefined(rc.params_, custom_params);
        default: return 0.0;
      }
    }

    // ========================================================================
    // CUDA Kernel
    // ========================================================================

    /// @brief Kernel: one thread per grid cell, iterates over all reactions
    __global__ void CalculateRateConstantsKernel(
        const double* __restrict__ d_temperature,
        const double* __restrict__ d_pressure,
        const double* __restrict__ d_air_density,
        const CudaMatrixParam custom_rate_params,
        const double* __restrict__ d_fixed_reactants,
        CudaMatrixParam rate_constants_param,
        const ProcessParam devstruct)
    {
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const std::size_t number_of_reactions = devstruct.number_of_reactions_;
      const std::size_t L = rate_constants_param.vector_length_;
      const std::size_t number_of_groups = (number_of_grid_cells + L - 1) / L;

      if (tid >= number_of_grid_cells)
        return;

      const std::size_t local_tid = tid % L;
      const std::size_t group_id = tid / L;

      double T = d_temperature[tid];
      double P = d_pressure[tid];
      double air_density = d_air_density[tid];

      // Pointer to rate_constants output for this group
      double* rc_base = rate_constants_param.d_data_ + group_id * (rate_constants_param.number_of_elements_ / number_of_groups);

      // Pointer to custom_rate_parameters for this group
      const double* custom_base = nullptr;
      if (custom_rate_params.d_data_ != nullptr)
        custom_base =
            custom_rate_params.d_data_ + group_id * (custom_rate_params.number_of_elements_ / number_of_groups);

      // Pointer to fixed_reactants for this group (same layout as rate_constants)
      const double* fixed_base = d_fixed_reactants + group_id * (rate_constants_param.number_of_elements_ / number_of_groups);

      std::size_t custom_offset = 0;
      for (std::size_t i_rxn = 0; i_rxn < number_of_reactions; ++i_rxn)
      {
        const CudaRateConstantData& rc = devstruct.rate_constants_[i_rxn];

        // Gather custom params for this cell (vectorized layout: col * L + local_tid)
        double custom_vals[2] = { 0.0, 0.0 };
        if (custom_base != nullptr && rc.number_of_custom_parameters_ > 0)
        {
          for (std::size_t i = 0; i < rc.number_of_custom_parameters_ && i < 2; ++i)
          {
            custom_vals[i] = custom_base[(custom_offset + i) * L + local_tid];
          }
        }

        double rate = CalculateRateConstant(rc, T, P, air_density, custom_vals);
        double fixed = fixed_base[i_rxn * L + local_tid];

        rc_base[i_rxn * L + local_tid] = rate * fixed;

        custom_offset += rc.number_of_custom_parameters_;
      }
    }

    // ========================================================================
    // Memory management functions
    // ========================================================================

    ProcessParam CopyProcessConstData(ProcessParam& hoststruct)
    {
      ProcessParam devstruct;
      devstruct.number_of_reactions_ = hoststruct.number_of_reactions_;

      std::size_t rate_constants_bytes = sizeof(CudaRateConstantData) * hoststruct.number_of_reactions_;

      auto cuda_stream_id = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.rate_constants_), rate_constants_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.rate_constants_,
              hoststruct.rate_constants_,
              rate_constants_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");

      return devstruct;
    }

    void FreeProcessConstData(ProcessParam& devstruct)
    {
      auto cuda_stream_id = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      if (devstruct.rate_constants_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.rate_constants_, cuda_stream_id), "cudaFree");
      devstruct.rate_constants_ = nullptr;
    }

    // ========================================================================
    // Kernel driver
    // ========================================================================

    void CalculateRateConstantsKernelDriver(
        const double* d_temperature,
        const double* d_pressure,
        const double* d_air_density,
        const CudaMatrixParam& custom_rate_params,
        const double* d_fixed_reactants,
        CudaMatrixParam& rate_constants_param,
        const ProcessParam& devstruct)
    {
      const std::size_t number_of_blocks = (rate_constants_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;

      CalculateRateConstantsKernel<<<
          number_of_blocks,
          BLOCK_SIZE,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          d_temperature, d_pressure, d_air_density, custom_rate_params, d_fixed_reactants, rate_constants_param, devstruct);
    }

  }  // namespace cuda
}  // namespace micm
