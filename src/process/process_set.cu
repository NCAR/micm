// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <micm/cuda/util/cuda_param.hpp>
#include <micm/cuda/util/cuda_util.cuh>

namespace micm
{
  namespace cuda
  {
    /// This is the CUDA kernel that calculates the forcing terms on the device
    __global__ void AddForcingTermsKernel(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam forcing_param,
        const ProcessSetParam devstruct)
    {
      // Calculate global thread ID
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      // Local device variables
      const std::size_t* const __restrict__ d_number_of_reactants = devstruct.number_of_reactants_;
      const std::size_t* __restrict__ d_reactant_ids = devstruct.reactant_ids_;
      const std::size_t* const __restrict__ d_number_of_products = devstruct.number_of_products_;
      const std::size_t* __restrict__ d_product_ids = devstruct.product_ids_;
      const double* __restrict__ d_yields = devstruct.yields_;
      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const double* __restrict__ d_rate_constants = rate_constants_param.d_data_;
      const double* __restrict__ d_state_variables = state_variables_param.d_data_;
      double* __restrict__ d_forcing = forcing_param.d_data_;

      const std::size_t cuda_matrix_vector_length = state_variables_param.vector_length_;
      const std::size_t number_of_groups =
          (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
      const std::size_t local_tid = tid % cuda_matrix_vector_length;
      const std::size_t group_id = tid / cuda_matrix_vector_length;
      const std::size_t number_of_reactions =
          rate_constants_param.number_of_elements_ / number_of_groups / cuda_matrix_vector_length;

      d_rate_constants += group_id * rate_constants_param.number_of_elements_ / number_of_groups;
      d_state_variables += group_id * state_variables_param.number_of_elements_ / number_of_groups;
      d_forcing += group_id * forcing_param.number_of_elements_ / number_of_groups;

      if (tid < number_of_grid_cells)
      {
        for (int i_rxn = 0; i_rxn < number_of_reactions; ++i_rxn)
        {
          double rate = d_rate_constants[(i_rxn * cuda_matrix_vector_length) + local_tid];
          const std::size_t number_of_reactants = d_number_of_reactants[i_rxn];
          for (int i_react = 0; i_react < number_of_reactants; ++i_react)
          {
            rate *= d_state_variables[(d_reactant_ids[i_react] * cuda_matrix_vector_length) + local_tid];
          }
          for (int i_react = 0; i_react < number_of_reactants; ++i_react)
          {
            d_forcing[(d_reactant_ids[i_react] * cuda_matrix_vector_length) + local_tid] -= rate;
          }
          const std::size_t number_of_products = d_number_of_products[i_rxn];
          for (int i_prod = 0; i_prod < number_of_products; ++i_prod)
          {
            d_forcing[d_product_ids[i_prod] * cuda_matrix_vector_length + local_tid] += d_yields[i_prod] * rate;
          }
          d_reactant_ids += d_number_of_reactants[i_rxn];
          d_product_ids += d_number_of_products[i_rxn];
          d_yields += d_number_of_products[i_rxn];
        }  // end of loop over number of reactions
      }  // end of checking a valid CUDA thread id
    }  // end of AddForcingTerms_kernel

    /// This is the CUDA kernel that forms the negative Jacobian matrix (-J) on the device
    __global__ void SubtractJacobianTermsKernel(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct)
    {
      // Calculate global thread ID
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      /// Local device variables
      const ProcessInfoParam* const __restrict__ d_jacobian_process_info = devstruct.jacobian_process_info_;
      const std::size_t* __restrict__ d_reactant_ids = devstruct.jacobian_reactant_ids_;
      const double* __restrict__ d_yields = devstruct.jacobian_yields_;
      const std::size_t* __restrict__ d_jacobian_flat_ids = devstruct.jacobian_flat_ids_;
      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const std::size_t number_of_process_infos = devstruct.jacobian_process_info_size_;
      const double* __restrict__ d_rate_constants = rate_constants_param.d_data_;
      const double* __restrict__ d_state_variables = state_variables_param.d_data_;
      double* __restrict__ d_jacobian = jacobian_param.d_data_;

      const std::size_t cuda_matrix_vector_length = state_variables_param.vector_length_;
      const std::size_t number_of_groups =
          (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
      const std::size_t local_tid = tid % cuda_matrix_vector_length;
      const std::size_t group_id = tid / cuda_matrix_vector_length;

      d_rate_constants += group_id * rate_constants_param.number_of_elements_ / number_of_groups;
      d_state_variables += group_id * state_variables_param.number_of_elements_ / number_of_groups;
      d_jacobian += group_id * jacobian_param.number_of_elements_ / number_of_groups;

      if (tid < number_of_grid_cells)
      {
        // loop over reactions in a grid
        for (int i_proc = 0; i_proc < number_of_process_infos; ++i_proc)
        {
          const ProcessInfoParam& process_info = d_jacobian_process_info[i_proc];
          // Calculate d_rate/d_ind
          double d_rate_d_ind = d_rate_constants[(process_info.process_id_ * cuda_matrix_vector_length) + local_tid];
          for (int i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
          {
            d_rate_d_ind *= d_state_variables[(d_reactant_ids[i_react] * cuda_matrix_vector_length) + local_tid];
          }
          for (int i_dep = 0; i_dep < process_info.number_of_dependent_reactants_ + 1; ++i_dep)
          {
            d_jacobian[*d_jacobian_flat_ids + local_tid] += d_rate_d_ind;
            ++d_jacobian_flat_ids;
          }
          for (int i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
          {
            std::size_t jacobian_idx = *d_jacobian_flat_ids + local_tid;
            d_jacobian[jacobian_idx] -= d_yields[i_dep] * d_rate_d_ind;
            ++d_jacobian_flat_ids;
          }
          d_reactant_ids += process_info.number_of_dependent_reactants_;
          d_yields += process_info.number_of_products_;
        }  // end of loop over reactions in a grid cell
      }  // end of checking a CUDA thread id
    }  // end of SubtractJacobianTermsKernel

    /// Unrolled variant of SubtractJacobianTermsKernel. Same data layout and access pattern
    /// as the original (pointer advancement over the const arrays); only difference is
    /// #pragma unroll hints on every loop.
    __global__ void SubtractJacobianTermsKernelUnrolled(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct)
    {
      std::size_t tid = blockIdx.x * BLOCK_SIZE + threadIdx.x;

      const ProcessInfoParam* const __restrict__ d_jacobian_process_info = devstruct.jacobian_process_info_;
      const std::size_t* __restrict__ d_reactant_ids = devstruct.jacobian_reactant_ids_;
      const double* __restrict__ d_yields = devstruct.jacobian_yields_;
      const std::size_t* __restrict__ d_jacobian_flat_ids = devstruct.jacobian_flat_ids_;
      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const std::size_t number_of_process_infos = devstruct.jacobian_process_info_size_;
      const double* __restrict__ d_rate_constants = rate_constants_param.d_data_;
      const double* __restrict__ d_state_variables = state_variables_param.d_data_;
      double* __restrict__ d_jacobian = jacobian_param.d_data_;

      const std::size_t cuda_matrix_vector_length = state_variables_param.vector_length_;
      const std::size_t number_of_groups =
          (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
      const std::size_t local_tid = tid % cuda_matrix_vector_length;
      const std::size_t group_id = tid / cuda_matrix_vector_length;

      d_rate_constants += group_id * rate_constants_param.number_of_elements_ / number_of_groups;
      d_state_variables += group_id * state_variables_param.number_of_elements_ / number_of_groups;
      d_jacobian += group_id * jacobian_param.number_of_elements_ / number_of_groups;

      if (tid < number_of_grid_cells)
      {
        for (int i_proc = 0; i_proc < number_of_process_infos; ++i_proc)
        {
          const ProcessInfoParam& process_info = d_jacobian_process_info[i_proc];
          double d_rate_d_ind = d_rate_constants[(process_info.process_id_ * cuda_matrix_vector_length) + local_tid];
          #pragma unroll
          for (int i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
          {
            d_rate_d_ind *= d_state_variables[(d_reactant_ids[i_react] * cuda_matrix_vector_length) + local_tid];
          }
          #pragma unroll
          for (int i_dep = 0; i_dep < process_info.number_of_dependent_reactants_ + 1; ++i_dep)
          {
            d_jacobian[*d_jacobian_flat_ids + local_tid] += d_rate_d_ind;
            ++d_jacobian_flat_ids;
          }
          #pragma unroll
          for (int i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
          {
            std::size_t jacobian_idx = *d_jacobian_flat_ids + local_tid;
            d_jacobian[jacobian_idx] -= d_yields[i_dep] * d_rate_d_ind;
            ++d_jacobian_flat_ids;
          }
          d_reactant_ids += process_info.number_of_dependent_reactants_;
          d_yields += process_info.number_of_products_;
        }
      }
    }  // end of SubtractJacobianTermsKernelUnrolled

    /// Shared-memory variant that parallelizes the i_proc loop across threads of a block.
    /// Layout: 1 block = 1 grid cell. blockDim.x threads cooperatively process the
    /// number_of_process_infos i_procs (each thread handles a contiguous chunk). Each
    /// thread accumulates its scatter contributions into a per-cell shared-memory
    /// accumulator using shared-memory atomicAdd (NO atomics on global d_jacobian).
    /// After a barrier, threads cooperatively write the accumulator out to global
    /// d_jacobian.
    __global__ void SubtractJacobianTermsKernelShared(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct)
    {
      extern __shared__ double s_jacobian[];

      const std::size_t cell_id = blockIdx.x;
      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      if (cell_id >= number_of_grid_cells)
        return;

      const std::size_t cuda_matrix_vector_length = state_variables_param.vector_length_;
      const std::size_t number_of_groups =
          (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
      const std::size_t group_id = cell_id / cuda_matrix_vector_length;
      const std::size_t local_tid = cell_id % cuda_matrix_vector_length;

      const std::size_t group_size_doubles = jacobian_param.number_of_elements_ / number_of_groups;
      const std::size_t flat_block_size = group_size_doubles / cuda_matrix_vector_length;

      const ProcessInfoParam* const __restrict__ d_jacobian_process_info = devstruct.jacobian_process_info_;
      const std::size_t* const __restrict__ d_reactant_ids_base = devstruct.jacobian_reactant_ids_;
      const double* const __restrict__ d_yields_base = devstruct.jacobian_yields_;
      const std::size_t* const __restrict__ d_flat_ids_base = devstruct.jacobian_flat_ids_;
      const std::size_t* const __restrict__ d_reactant_ids_offsets = devstruct.jacobian_reactant_ids_offsets_;
      const std::size_t* const __restrict__ d_yields_offsets = devstruct.jacobian_yields_offsets_;
      const std::size_t* const __restrict__ d_flat_ids_offsets = devstruct.jacobian_flat_ids_offsets_;
      const std::size_t number_of_process_infos = devstruct.jacobian_process_info_size_;

      const double* __restrict__ d_rate_constants =
          rate_constants_param.d_data_ + group_id * (rate_constants_param.number_of_elements_ / number_of_groups);
      const double* __restrict__ d_state_variables =
          state_variables_param.d_data_ + group_id * (state_variables_param.number_of_elements_ / number_of_groups);
      double* __restrict__ d_jacobian_group = jacobian_param.d_data_ + group_id * group_size_doubles;

      // Phase 1: zero shared accumulator
      for (std::size_t i = threadIdx.x; i < flat_block_size; i += blockDim.x)
        s_jacobian[i] = 0.0;
      __syncthreads();

      // Phase 2: parallelize i_proc loop. Each thread takes a contiguous chunk so that the
      // precomputed offsets give it a direct entry point into the const slices.
      const std::size_t chunk = (number_of_process_infos + blockDim.x - 1) / blockDim.x;
      const std::size_t my_start = static_cast<std::size_t>(threadIdx.x) * chunk;
      std::size_t my_end = my_start + chunk;
      if (my_end > number_of_process_infos)
        my_end = number_of_process_infos;

      for (std::size_t i_proc = my_start; i_proc < my_end; ++i_proc)
      {
        const ProcessInfoParam process_info = d_jacobian_process_info[i_proc];
        const std::size_t* my_reactant_ids = d_reactant_ids_base + d_reactant_ids_offsets[i_proc];
        const double* my_yields = d_yields_base + d_yields_offsets[i_proc];
        const std::size_t* my_flat_ids = d_flat_ids_base + d_flat_ids_offsets[i_proc];

        double d_rate_d_ind = d_rate_constants[(process_info.process_id_ * cuda_matrix_vector_length) + local_tid];
        for (std::size_t i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
        {
          d_rate_d_ind *= d_state_variables[(my_reactant_ids[i_react] * cuda_matrix_vector_length) + local_tid];
        }

        // Scatter into shared accumulator. Multiple threads of this block may target the
        // same flat_id, so use shared-memory atomicAdd. This is NOT a global-memory atomic.
        // jacobian_flat_ids values are pre-scaled by cuda_matrix_vector_length (they index
        // into the per-element group buffer); divide them out to get the per-row index
        // into s_jacobian (which holds one slot per sparse entry, not per element).
        for (std::size_t i_dep = 0; i_dep < process_info.number_of_dependent_reactants_ + 1; ++i_dep)
        {
          atomicAdd(&s_jacobian[my_flat_ids[i_dep] / cuda_matrix_vector_length], d_rate_d_ind);
        }
        const std::size_t prod_offset = process_info.number_of_dependent_reactants_ + 1;
        for (std::size_t i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
        {
          atomicAdd(
              &s_jacobian[my_flat_ids[prod_offset + i_dep] / cuda_matrix_vector_length],
              -my_yields[i_dep] * d_rate_d_ind);
        }
      }
      __syncthreads();

      // Phase 3: accumulate this cell's slice into d_jacobian. Use += to match the
      // original kernel, which accumulates across invocations (the caller may issue
      // multiple launches between memsets). Different blocks within the same group
      // target distinct local_tid columns, so no inter-block race here.
      // (uncoalesced: stride = cuda_matrix_vector_length within each warp; see plan trade-off).
      for (std::size_t i = threadIdx.x; i < flat_block_size; i += blockDim.x)
      {
        d_jacobian_group[i * cuda_matrix_vector_length + local_tid] += s_jacobian[i];
      }
    }  // end of SubtractJacobianTermsKernelShared

    /// Strided-i_proc, atomic-free variant. Each thread strides across the i_proc loop
    /// (i_proc = threadIdx.x, threadIdx.x + blockDim.x, ...), accumulates its
    /// (slot, value) contributions into a per-thread register-resident list, and the
    /// block then merges those lists into a per-cell s_jacobian via a chunked
    /// shared-memory tree reduction (warp-shuffle butterfly within each slot column).
    /// No atomics anywhere — neither global nor shared. Reduction order is fixed
    /// (thread 0 -> blockDim.x-1 ascending), so the result is deterministic across runs.
    __global__ void SubtractJacobianTermsKernelSharedReduce(
        const CudaMatrixParam rate_constants_param,
        const CudaMatrixParam state_variables_param,
        CudaMatrixParam jacobian_param,
        const ProcessSetParam devstruct)
    {
      // Per-thread upper bound on contributions: ceil(JACOBIAN_FLAT_IDS_SIZE / BLOCK_SIZE).
      // For TS1 that's ceil(4255/128) = 34; pad to 48 for safety/headroom.
      constexpr int MAX_CONTRIB = 48;
      // Slots reduced per chunk pass. Smem cost = CHUNK_SIZE * BLOCK_SIZE * sizeof(double).
      // 16 * 128 * 8 = 16,384 B; plus s_jacobian (FLAT_BLOCK_SIZE doubles) we stay <48 KB.
      constexpr int CHUNK_SIZE = 16;

      extern __shared__ double s_mem[];

      const std::size_t cell_id = blockIdx.x;
      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      if (cell_id >= number_of_grid_cells)
        return;

      const std::size_t cuda_matrix_vector_length = state_variables_param.vector_length_;
      const std::size_t number_of_groups =
          (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
      const std::size_t group_id = cell_id / cuda_matrix_vector_length;
      const std::size_t local_tid = cell_id % cuda_matrix_vector_length;

      const std::size_t group_size_doubles = jacobian_param.number_of_elements_ / number_of_groups;
      const std::size_t flat_block_size = group_size_doubles / cuda_matrix_vector_length;

      const ProcessInfoParam* const __restrict__ d_jacobian_process_info = devstruct.jacobian_process_info_;
      const std::size_t* const __restrict__ d_reactant_ids_base = devstruct.jacobian_reactant_ids_;
      const double* const __restrict__ d_yields_base = devstruct.jacobian_yields_;
      const std::size_t* const __restrict__ d_flat_ids_base = devstruct.jacobian_flat_ids_;
      const std::size_t* const __restrict__ d_reactant_ids_offsets = devstruct.jacobian_reactant_ids_offsets_;
      const std::size_t* const __restrict__ d_yields_offsets = devstruct.jacobian_yields_offsets_;
      const std::size_t* const __restrict__ d_flat_ids_offsets = devstruct.jacobian_flat_ids_offsets_;
      const std::size_t number_of_process_infos = devstruct.jacobian_process_info_size_;

      const double* __restrict__ d_rate_constants =
          rate_constants_param.d_data_ + group_id * (rate_constants_param.number_of_elements_ / number_of_groups);
      const double* __restrict__ d_state_variables =
          state_variables_param.d_data_ + group_id * (state_variables_param.number_of_elements_ / number_of_groups);
      double* __restrict__ d_jacobian_group = jacobian_param.d_data_ + group_id * group_size_doubles;

      // Smem layout: [s_jacobian (flat_block_size doubles)] [s_buf (CHUNK_SIZE * blockDim.x doubles)]
      double* __restrict__ s_jacobian = s_mem;
      double* __restrict__ s_buf = s_mem + flat_block_size;

      // Phase 0: zero s_jacobian.
      for (std::size_t i = threadIdx.x; i < flat_block_size; i += blockDim.x)
        s_jacobian[i] = 0.0;
      __syncthreads();

      // Phase 1: strided i_proc walk; build per-thread (slot, value) list in registers.
      // Slots fit in 16 bits since flat_block_size < 65536 for all supported mechanisms.
      uint16_t my_slot[MAX_CONTRIB];
      double my_value[MAX_CONTRIB];
      int my_n = 0;

      for (std::size_t i_proc = threadIdx.x; i_proc < number_of_process_infos; i_proc += blockDim.x)
      {
        const ProcessInfoParam process_info = d_jacobian_process_info[i_proc];
        const std::size_t* my_reactant_ids = d_reactant_ids_base + d_reactant_ids_offsets[i_proc];
        const double* my_yields = d_yields_base + d_yields_offsets[i_proc];
        const std::size_t* my_flat_ids = d_flat_ids_base + d_flat_ids_offsets[i_proc];

        double d_rate_d_ind = d_rate_constants[(process_info.process_id_ * cuda_matrix_vector_length) + local_tid];
        for (std::size_t i_react = 0; i_react < process_info.number_of_dependent_reactants_; ++i_react)
        {
          d_rate_d_ind *= d_state_variables[(my_reactant_ids[i_react] * cuda_matrix_vector_length) + local_tid];
        }

        const std::size_t n_pos = process_info.number_of_dependent_reactants_ + 1;
        for (std::size_t i_dep = 0; i_dep < n_pos; ++i_dep)
        {
          my_slot[my_n] = static_cast<uint16_t>(my_flat_ids[i_dep] / cuda_matrix_vector_length);
          my_value[my_n] = d_rate_d_ind;
          ++my_n;
        }
        for (std::size_t i_dep = 0; i_dep < process_info.number_of_products_; ++i_dep)
        {
          my_slot[my_n] = static_cast<uint16_t>(my_flat_ids[n_pos + i_dep] / cuda_matrix_vector_length);
          my_value[my_n] = -my_yields[i_dep] * d_rate_d_ind;
          ++my_n;
        }
      }

      // Pre-sort per-thread list by slot (insertion sort, O(my_n^2) ≈ 1k ops worst case).
      // After sorting, Phase 2b can walk linearly with a monotonic cursor instead of
      // scanning all my_n entries per chunk.
      for (int i = 1; i < my_n; ++i)
      {
        const uint16_t s = my_slot[i];
        const double v = my_value[i];
        int j = i - 1;
        while (j >= 0 && my_slot[j] > s)
        {
          my_slot[j + 1] = my_slot[j];
          my_value[j + 1] = my_value[j];
          --j;
        }
        my_slot[j + 1] = s;
        my_value[j + 1] = v;
      }

      // Phase 2: chunked atomic-free reduction into s_jacobian.
      const int warp_id = threadIdx.x / 32;
      const int lane = threadIdx.x % 32;
      const int warps_per_block = blockDim.x / 32;
      int cursor = 0;  // monotonic into the sorted per-thread list

      for (std::size_t slot_base = 0; slot_base < flat_block_size; slot_base += CHUNK_SIZE)
      {
        // 2a) Zero this chunk's portion of s_buf.
        for (std::size_t k = threadIdx.x; k < CHUNK_SIZE * blockDim.x; k += blockDim.x)
          s_buf[k] = 0.0;
        __syncthreads();

        // 2b) Drop this thread's contributions whose slot lies in [slot_base, slot_base+CHUNK_SIZE)
        //     into s_buf[slot_in_chunk * blockDim.x + threadIdx.x]. Single-writer per cell — no race.
        //     Linear walk thanks to the sort above.
        const std::size_t slot_end = slot_base + CHUNK_SIZE;
        while (cursor < my_n && static_cast<std::size_t>(my_slot[cursor]) < slot_end)
        {
          const int s = static_cast<int>(my_slot[cursor]) - static_cast<int>(slot_base);
          s_buf[s * blockDim.x + threadIdx.x] += my_value[cursor];
          ++cursor;
        }
        __syncthreads();

        // 2c) Reduce CHUNK_SIZE slot columns of length blockDim.x. One warp per slot
        //     (warps_per_block warps handle CHUNK_SIZE slots, striding by warps_per_block).
        //     Within a warp: each lane sums blockDim.x/32 strided values in ascending order,
        //     then a 5-step __shfl_down_sync butterfly. Lane 0 writes the single += into s_jacobian.
        for (int k = warp_id; k < CHUNK_SIZE; k += warps_per_block)
        {
          double v = 0.0;
          // For BLOCK_SIZE=128, this is 4 strided loads in ascending order: lane, lane+32, lane+64, lane+96.
          for (int t = lane; t < static_cast<int>(blockDim.x); t += 32)
            v += s_buf[k * blockDim.x + t];

          v += __shfl_down_sync(0xffffffff, v, 16);
          v += __shfl_down_sync(0xffffffff, v, 8);
          v += __shfl_down_sync(0xffffffff, v, 4);
          v += __shfl_down_sync(0xffffffff, v, 2);
          v += __shfl_down_sync(0xffffffff, v, 1);

          if (lane == 0 && (slot_base + k) < flat_block_size)
            s_jacobian[slot_base + k] += v;
        }
        __syncthreads();
      }

      // Phase 3: drain s_jacobian into the per-cell column of d_jacobian with += semantics
      // to match the original kernel (which accumulates across invocations between memsets).
      for (std::size_t i = threadIdx.x; i < flat_block_size; i += blockDim.x)
      {
        d_jacobian_group[i * cuda_matrix_vector_length + local_tid] += s_jacobian[i];
      }
    }  // end of SubtractJacobianTermsKernelSharedReduce

    /// This is the function that will copy the constant data
    ///   members of class "CudaProcessSet" to the device,
    ///   except for the "jacobian_flat_id" because it is unknown now
    ProcessSetParam CopyConstData(ProcessSetParam& hoststruct)
    {
      /// Calculate the memory space of each constant data member
      std::size_t number_of_reactants_bytes = sizeof(std::size_t) * hoststruct.number_of_reactants_size_;
      std::size_t reactant_ids_bytes = sizeof(std::size_t) * hoststruct.reactant_ids_size_;
      std::size_t number_of_products_bytes = sizeof(std::size_t) * hoststruct.number_of_products_size_;
      std::size_t product_ids_bytes = sizeof(std::size_t) * hoststruct.product_ids_size_;
      std::size_t yields_bytes = sizeof(double) * hoststruct.yields_size_;

      /// Create a struct whose members contain the addresses in the device memory.
      ProcessSetParam devstruct;

      auto cuda_stream_id = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      /// Allocate memory space on the device
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.number_of_reactants_), number_of_reactants_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.reactant_ids_), reactant_ids_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.number_of_products_), number_of_products_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.product_ids_), product_ids_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.yields_), yields_bytes, cuda_stream_id), "cudaMalloc");

      /// Copy the data from host to device
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.number_of_reactants_,
              hoststruct.number_of_reactants_,
              number_of_reactants_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.reactant_ids_, hoststruct.reactant_ids_, reactant_ids_bytes, cudaMemcpyHostToDevice, cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.number_of_products_,
              hoststruct.number_of_products_,
              number_of_products_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.product_ids_, hoststruct.product_ids_, product_ids_bytes, cudaMemcpyHostToDevice, cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(devstruct.yields_, hoststruct.yields_, yields_bytes, cudaMemcpyHostToDevice, cuda_stream_id),
          "cudaMemcpy");

      devstruct.number_of_reactants_size_ = hoststruct.number_of_reactants_size_;
      devstruct.reactant_ids_size_ = hoststruct.reactant_ids_size_;
      devstruct.number_of_products_size_ = hoststruct.number_of_products_size_;
      devstruct.product_ids_size_ = hoststruct.product_ids_size_;
      devstruct.yields_size_ = hoststruct.yields_size_;

      return devstruct;
    }

    /// This is the function that will copy the information needed
    /// to calculate Jacobian terms to the device
    void CopyJacobianParams(ProcessSetParam& hoststruct, ProcessSetParam& devstruct)
    {
      /// Calculate the memory space
      std::size_t jacobian_process_info_bytes = sizeof(ProcessInfoParam) * hoststruct.jacobian_process_info_size_;
      std::size_t jacobian_reactant_ids_bytes = sizeof(std::size_t) * hoststruct.jacobian_reactant_ids_size_;
      std::size_t jacobian_product_ids_bytes = sizeof(std::size_t) * hoststruct.jacobian_product_ids_size_;
      std::size_t jacobian_yields_bytes = sizeof(double) * hoststruct.jacobian_yields_size_;
      std::size_t jacobian_flat_ids_bytes = sizeof(std::size_t) * hoststruct.jacobian_flat_ids_size_;

      auto cuda_stream_id = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      /// Allocate memory space on the device
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jacobian_process_info_), jacobian_process_info_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jacobian_reactant_ids_), jacobian_reactant_ids_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jacobian_product_ids_), jacobian_product_ids_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(cudaMallocAsync(&(devstruct.jacobian_yields_), jacobian_yields_bytes, cuda_stream_id), "cudaMalloc");
      CHECK_CUDA_ERROR(
          cudaMallocAsync(&(devstruct.jacobian_flat_ids_), jacobian_flat_ids_bytes, cuda_stream_id), "cudaMalloc");

      /// Copy the data from host to device
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_process_info_,
              hoststruct.jacobian_process_info_,
              jacobian_process_info_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_reactant_ids_,
              hoststruct.jacobian_reactant_ids_,
              jacobian_reactant_ids_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_product_ids_,
              hoststruct.jacobian_product_ids_,
              jacobian_product_ids_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_yields_,
              hoststruct.jacobian_yields_,
              jacobian_yields_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");
      CHECK_CUDA_ERROR(
          cudaMemcpyAsync(
              devstruct.jacobian_flat_ids_,
              hoststruct.jacobian_flat_ids_,
              jacobian_flat_ids_bytes,
              cudaMemcpyHostToDevice,
              cuda_stream_id),
          "cudaMemcpy");

      devstruct.jacobian_process_info_size_ = hoststruct.jacobian_process_info_size_;
      devstruct.jacobian_reactant_ids_size_ = hoststruct.jacobian_reactant_ids_size_;
      devstruct.jacobian_product_ids_size_ = hoststruct.jacobian_product_ids_size_;
      devstruct.jacobian_yields_size_ = hoststruct.jacobian_yields_size_;
      devstruct.jacobian_flat_ids_size_ = hoststruct.jacobian_flat_ids_size_;
    }

    /// This is the function that will delete the constant data
    ///   members of class "CudaProcessSet" on the device
    void FreeConstData(ProcessSetParam& devstruct)
    {
      auto cuda_stream_id = micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0);

      if (devstruct.number_of_reactants_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.number_of_reactants_, cuda_stream_id), "cudaFree");
      if (devstruct.reactant_ids_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.reactant_ids_, cuda_stream_id), "cudaFree");
      if (devstruct.number_of_products_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.number_of_products_, cuda_stream_id), "cudaFree");
      if (devstruct.product_ids_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.product_ids_, cuda_stream_id), "cudaFree");
      if (devstruct.yields_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.yields_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_process_info_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_process_info_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_reactant_ids_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_reactant_ids_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_product_ids_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_product_ids_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_yields_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_yields_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_flat_ids_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_flat_ids_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_reactant_ids_offsets_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_reactant_ids_offsets_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_yields_offsets_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_yields_offsets_, cuda_stream_id), "cudaFree");
      if (devstruct.jacobian_flat_ids_offsets_ != nullptr)
        CHECK_CUDA_ERROR(cudaFreeAsync(devstruct.jacobian_flat_ids_offsets_, cuda_stream_id), "cudaFree");
    }

    void SubtractJacobianTermsKernelDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      const std::size_t number_of_blocks = (rate_constants_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      SubtractJacobianTermsKernel<<<
          number_of_blocks,
          BLOCK_SIZE,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, jacobian_param, devstruct);
    }  // end of SubtractJacobianTermsKernelDriver

    void AddForcingTermsKernelDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& forcing_param,
        const ProcessSetParam& devstruct)
    {
      const std::size_t number_of_blocks = (rate_constants_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      AddForcingTermsKernel<<<
          number_of_blocks,
          BLOCK_SIZE,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, forcing_param, devstruct);
    }  // end of AddForcingTermsKernelDriver

    void SubtractJacobianTermsKernelUnrolledDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      const std::size_t number_of_blocks = (rate_constants_param.number_of_grid_cells_ + BLOCK_SIZE - 1) / BLOCK_SIZE;
      SubtractJacobianTermsKernelUnrolled<<<
          number_of_blocks,
          BLOCK_SIZE,
          0,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, jacobian_param, devstruct);
    }  // end of SubtractJacobianTermsKernelUnrolledDriver

    void SubtractJacobianTermsKernelSharedDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const std::size_t cuda_matrix_vector_length = state_variables_param.vector_length_;
      const std::size_t number_of_groups =
          (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
      const std::size_t flat_block_size =
          jacobian_param.number_of_elements_ / number_of_groups / cuda_matrix_vector_length;
      const std::size_t smem_bytes = flat_block_size * sizeof(double);

      SubtractJacobianTermsKernelShared<<<
          number_of_grid_cells,
          BLOCK_SIZE,
          smem_bytes,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, jacobian_param, devstruct);
    }  // end of SubtractJacobianTermsKernelSharedDriver

    void SubtractJacobianTermsKernelSharedReduceDriver(
        const CudaMatrixParam& rate_constants_param,
        const CudaMatrixParam& state_variables_param,
        CudaMatrixParam& jacobian_param,
        const ProcessSetParam& devstruct)
    {
      const std::size_t number_of_grid_cells = rate_constants_param.number_of_grid_cells_;
      const std::size_t cuda_matrix_vector_length = state_variables_param.vector_length_;
      const std::size_t number_of_groups =
          (number_of_grid_cells + cuda_matrix_vector_length - 1) / cuda_matrix_vector_length;
      const std::size_t flat_block_size =
          jacobian_param.number_of_elements_ / number_of_groups / cuda_matrix_vector_length;

      // Mirror of the kernel-side CHUNK_SIZE constant for sizing the smem allocation.
      constexpr std::size_t CHUNK_SIZE = 16;
      const std::size_t smem_bytes = sizeof(double) * (flat_block_size + CHUNK_SIZE * BLOCK_SIZE);

      SubtractJacobianTermsKernelSharedReduce<<<
          number_of_grid_cells,
          BLOCK_SIZE,
          smem_bytes,
          micm::cuda::CudaStreamSingleton::GetInstance().GetCudaStream(0)>>>(
          rate_constants_param, state_variables_param, jacobian_param, devstruct);
    }  // end of SubtractJacobianTermsKernelSharedReduceDriver
  }  // namespace cuda
}  // namespace micm
