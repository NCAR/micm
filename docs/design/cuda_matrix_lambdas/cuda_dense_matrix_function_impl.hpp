// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_dense_matrix_function.hpp>

namespace micm
{
  namespace cuda
  {
    // ============================================================================
    // CUDA Dense Matrix Function() Implementation Sketch
    // ============================================================================

    /**
     * @brief Callable wrapper that launches compiled CUDA kernel
     * 
     * This is returned by CudaDenseMatrix::Function() and can be invoked
     * multiple times with different matrices (of the same dimensions).
     */
    template<typename... MatrixTypes>
    class CompiledCudaFunction
    {
     private:
      CUfunction kernel_;
      OperationSequence ops_;
      std::vector<std::pair<std::size_t, std::size_t>> expected_dims_;  // (rows, cols) for each matrix
      std::size_t vector_length_;
      
     public:
      CompiledCudaFunction(
          CUfunction kernel,
          const OperationSequence& ops,
          const std::vector<std::pair<std::size_t, std::size_t>>& dims,
          std::size_t vec_len)
          : kernel_(kernel), ops_(ops), expected_dims_(dims), vector_length_(vec_len)
      {
      }
      
      /// @brief Execute the compiled kernel on the provided matrices
      void operator()(MatrixTypes&... matrices)
      {
        // Validate dimensions match expected
        std::size_t idx = 0;
        ([&](auto& matrix) {
          if (matrix.NumRows() != expected_dims_[idx].first ||
              matrix.NumColumns() != expected_dims_[idx].second)
          {
            throw std::system_error(
                make_error_code(MicmMatrixErrc::InvalidVector),
                "Matrix dimensions do not match compiled kernel");
          }
          ++idx;
        }(matrices), ...);
        
        // Ensure data is on device
        (matrices.CopyToDevice(), ...);
        
        // Set up kernel parameters
        std::vector<void*> kernel_params;
        ([&](auto& matrix) {
          auto param = matrix.AsDeviceParam();
          kernel_params.push_back(&param.d_data_);
          std::size_t rows = matrix.NumRows();
          std::size_t cols = matrix.NumColumns();
          kernel_params.push_back(&rows);
          kernel_params.push_back(&cols);
        }(matrices), ...);
        
        kernel_params.push_back(&vector_length_);
        
        // Calculate grid/block dimensions
        std::size_t total_elements = expected_dims_[0].first * vector_length_;
        int blockSize = 256;
        int numBlocks = (total_elements + blockSize - 1) / blockSize;
        
        // Launch kernel
        CUresult result = cuLaunchKernel(
            kernel_,
            numBlocks, 1, 1,    // grid dim
            blockSize, 1, 1,    // block dim
            0,                   // shared memory
            nullptr,             // stream
            kernel_params.data(),
            nullptr);
        
        if (result != CUDA_SUCCESS)
        {
          throw std::runtime_error("CUDA kernel launch failed");
        }
        
        // Wait for completion
        cuStreamSynchronize(nullptr);
        
        // Note: User must call CopyToHost() explicitly if they want results back
      }
    };

    // ============================================================================
    // Template Magic: Converting Lambda Calls to Recordings
    // ============================================================================

    /**
     * @brief Helper to extract lambda parameter types and execute with recording proxies
     * 
     * This is the trickiest part - we need to intercept the lambda's execution
     * and replace actual values with our recording proxies.
     */
    template<typename Func, typename... Args>
    class LambdaRecorder
    {
     public:
      /**
       * The key insight: When ForEachRow is called with a lambda like:
       *   [&](const double& a, const double& b, double& t) { t = a + b; }
       * 
       * We want to execute it with RecordingValue proxies instead of doubles:
       *   [&](RecordingValue& a, RecordingValue& b, RecordingValue& t) { t = a + b; }
       * 
       * The assignment `t = a + b` will record the expression tree.
       */
      
      static OperationSequence::ForEachRowCall record(
          RecordingContext* ctx,
          Func&& func,
          Args&&... args)
      {
        OperationSequence::ForEachRowCall call;
        
        // Convert each argument to a RecordingValue
        auto recording_args = convertArgsToRecording(ctx, std::forward<Args>(args)...);
        
        // Execute the lambda with recording arguments
        // The lambda body will execute, building the expression tree
        std::apply([&](auto&&... rec_args) {
          func(rec_args.asValue()...);
        }, recording_args);
        
        // Extract the recorded operations from context
        // TODO: Populate call.statements, call.input_columns, etc.
        
        return call;
      }
      
     private:
      // Convert ColumnView/RowVariable to recording proxies
      template<typename Arg>
      static auto convertToRecording(RecordingContext* ctx, Arg&& arg)
      {
        if constexpr (std::is_same_v<std::decay_t<Arg>, RecordingColumnView>)
        {
          return arg;  // Already a recording proxy
        }
        else if constexpr (std::is_same_v<std::decay_t<Arg>, RecordingRowVariable>)
        {
          return arg;  // Already a recording proxy
        }
        else
        {
          // Should not reach here in normal use
          static_assert(sizeof(Arg) == 0, "Unsupported argument type for recording");
        }
      }
      
      template<typename... ConvertArgs>
      static auto convertArgsToRecording(RecordingContext* ctx, ConvertArgs&&... args)
      {
        return std::make_tuple(convertToRecording(ctx, std::forward<ConvertArgs>(args))...);
      }
    };

  }  // namespace cuda

  // ============================================================================
  // CudaDenseMatrix::Function() Implementation
  // ============================================================================

  template<class T, std::size_t L>
  class CudaDenseMatrix;  // Forward declaration

  // Partial specialization or extension to add Function() method
  template<class T, std::size_t L>
  class CudaDenseMatrix
  {
   public:
    // ... existing methods ...
    
    /**
     * @brief Create a JIT-compiled CUDA kernel from a lambda function
     * 
     * This method records the operations specified in the lambda, generates
     * CUDA kernel code, compiles it using NVRTC, and returns a callable object
     * that launches the compiled kernel.
     * 
     * @param func Lambda function describing matrix operations
     * @param matrices Matrices to operate on (used for dimension validation)
     * @return Callable object that executes the compiled kernel
     * 
     * @example
     * ```cpp
     * CudaDenseMatrix<double, 4> A{3, 2};
     * CudaDenseMatrix<double, 4> B{3, 3};
     * 
     * auto func = CudaDenseMatrix<double, 4>::Function(
     *     [](auto&& mA, auto&& mB) {
     *         auto tmp = mA.GetRowVariable();
     *         mA.ForEachRow([&](const double& a, const double& b, double& t) {
     *             t = a + b;
     *         }, mA.GetConstColumnView(0), mB.GetConstColumnView(2), tmp);
     *         mA.ForEachRow([&](const double& t, double& c) {
     *             c = t;
     *         }, tmp, mA.GetColumnView(1));
     *     }, A, B);
     * 
     * func(A, B);  // Executes on GPU
     * ```
     */
    template<typename Func, typename... Matrices>
    static auto Function(Func&& func, Matrices&... matrices)
    {
      // Phase 1: Validate all matrices have same number of rows
      std::size_t num_rows = 0;
      std::vector<std::pair<std::size_t, std::size_t>> dims;
      std::size_t index = 0;
      
      ([&](auto& matrix) {
        if (index == 0)
        {
          num_rows = matrix.NumRows();
        }
        else if (matrix.NumRows() != num_rows)
        {
          throw std::system_error(
              make_error_code(MicmMatrixErrc::InvalidVector),
              "All matrices must have the same number of rows");
        }
        dims.push_back({matrix.NumRows(), matrix.NumColumns()});
        ++index;
      }(matrices), ...);
      
      // Phase 2: Create recording context and execute lambda with recording proxies
      cuda::RecordingContext ctx;
      ctx.sequence_.num_matrices = sizeof...(Matrices);
      ctx.sequence_.matrix_dimensions = dims;
      
      // Create recording group views for each matrix
      index = 0;
      auto recording_matrices = std::make_tuple(
          ([&](auto& /* matrix */) {
            return cuda::RecordingGroupView<L>(&ctx, index++);
          }(matrices))...);
      
      // Execute the user's lambda with recording proxies
      // This builds the operation sequence in ctx
      std::apply([&](auto&&... rec_mats) {
        func(rec_mats...);
      }, recording_matrices);
      
      // Phase 3: Generate CUDA kernel source code
      std::string kernel_source = cuda::CudaCodeGenerator::generateKernel(ctx.sequence_);
      
      // Optional: Print generated kernel for debugging
      #ifdef MICM_CUDA_FUNCTION_DEBUG
      std::cout << "Generated CUDA kernel:\n" << kernel_source << std::endl;
      #endif
      
      // Phase 4: Compile kernel (or retrieve from cache)
      CUfunction kernel = cuda::KernelCache::instance().getKernel(kernel_source);
      
      // Phase 5: Return callable wrapper
      return cuda::CompiledCudaFunction<Matrices...>(kernel, ctx.sequence_, dims, L);
    }
  };

}  // namespace micm
