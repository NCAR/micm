// Example: Using Function() with CPU and CUDA Matrices
// This demonstrates the intended API compatibility between CPU and GPU

#include <micm/util/vector_matrix.hpp>
#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <micm/cuda/util/cuda_dense_matrix_function.hpp>

#include <iostream>

int main()
{
  constexpr size_t L = 4;  // Vector size
  
  // ============================================================================
  // CPU VERSION - This works today
  // ============================================================================
  
  {
    using Matrix = micm::VectorMatrix<double, L>;
    
    Matrix matrixA{3, 2, 1.0};
    Matrix matrixB{3, 3, 2.0};
    
    // Initialize matrices
    for (int i = 0; i < 3; ++i)
    {
      matrixA[i][0] = static_cast<double>(i);
      matrixB[i][2] = static_cast<double>(i * 4);
    }
    
    // Create reusable function
    auto func = Matrix::Function(
        [](auto&& mA, auto&& mB)
        {
          // tmp = mA.col(0) + mB.col(2)
          auto tmp = mA.GetRowVariable();
          mA.ForEachRow(
              [&](const double& a, const double& b, double& t) { t = a + b; },
              mA.GetConstColumnView(0),
              mB.GetConstColumnView(2),
              tmp);
          
          // mA.col(1) = tmp
          mA.ForEachRow(
              [&](const double& t, double& c) { c = t; },
              tmp,
              mA.GetColumnView(1));
        },
        matrixA,
        matrixB);
    
    // Execute function
    func(matrixA, matrixB);
    
    // Check results
    std::cout << "CPU Results:\n";
    for (int i = 0; i < 3; ++i)
    {
      std::cout << "  Row " << i << ": " << matrixA[i][1] 
                << " (expected: " << (i + i*4) << ")\n";
    }
  }
  
  // ============================================================================
  // CUDA VERSION - Proposed implementation
  // ============================================================================
  
  {
    using Matrix = micm::CudaDenseMatrix<double, L>;
    
    Matrix matrixA{3, 2, 1.0};
    Matrix matrixB{3, 3, 2.0};
    
    // Initialize matrices on host
    for (int i = 0; i < 3; ++i)
    {
      matrixA[i][0] = static_cast<double>(i);
      matrixB[i][2] = static_cast<double>(i * 4);
    }
    
    // Copy to device
    matrixA.CopyToDevice();
    matrixB.CopyToDevice();
    
    // Create reusable function
    // NOTE: With proposed implementation, requires 'auto' parameters in nested lambdas
    auto func = Matrix::Function(
        [](auto&& mA, auto&& mB)
        {
          // tmp = mA.col(0) + mB.col(2)
          auto tmp = mA.GetRowVariable();
          mA.ForEachRow(
              [&](auto a, auto b, auto& t) { t = a + b; },  // 'auto' instead of 'const double&'
              mA.GetConstColumnView(0),
              mB.GetConstColumnView(2),
              tmp);
          
          // mA.col(1) = tmp
          mA.ForEachRow(
              [&](auto t, auto& c) { c = t; },  // 'auto' instead of 'const double&'
              tmp,
              mA.GetColumnView(1));
        },
        matrixA,  // Used for dimension inference
        matrixB);
    
    // On first call: Compiles CUDA kernel (~200-500ms)
    // Subsequent calls: Launches cached kernel (~microseconds)
    func(matrixA, matrixB);
    
    // Copy results back to host
    matrixA.CopyToHost();
    
    // Check results
    std::cout << "\nCUDA Results:\n";
    for (int i = 0; i < 3; ++i)
    {
      std::cout << "  Row " << i << ": " << matrixA[i][1] 
                << " (expected: " << (i + i*4) << ")\n";
    }
  }
  
  // ============================================================================
  // More Complex Example: Multiple Temporaries and Math Functions
  // ============================================================================
  
  {
    using Matrix = micm::CudaDenseMatrix<double, L>;
    
    Matrix matrix{4, 5, 1.0};
    
    // Initialize
    for (size_t i = 0; i < matrix.NumRows(); ++i)
    {
      matrix[i][0] = static_cast<double>(i + 1);
      matrix[i][1] = static_cast<double>((i + 1) * 10);
    }
    
    matrix.CopyToDevice();
    
    // Complex operations
    auto func = Matrix::Function(
        [](auto&& m)
        {
          auto tmp1 = m.GetRowVariable();
          auto tmp2 = m.GetRowVariable();
          
          // tmp1 = col0 * col1
          m.ForEachRow(
              [&](auto a, auto b, auto& t) { t = a * b; },
              m.GetConstColumnView(0),
              m.GetConstColumnView(1),
              tmp1);
          
          // tmp2 = exp(col0) + log(col1)
          m.ForEachRow(
              [&](auto a, auto b, auto& t) { 
                  t = exp(a) + log(b);
              },
              m.GetConstColumnView(0),
              m.GetConstColumnView(1),
              tmp2);
          
          // col2 = tmp1 + tmp2
          m.ForEachRow(
              [&](auto t1, auto t2, auto& c) { c = t1 + t2; },
              tmp1,
              tmp2,
              m.GetColumnView(2));
          
          // col3 = sqrt(col2) / col0
          m.ForEachRow(
              [&](auto c2, auto c0, auto& c3) { 
                  c3 = sqrt(c2) / c0; 
              },
              m.GetConstColumnView(2),
              m.GetConstColumnView(0),
              m.GetColumnView(3));
        },
        matrix);
    
    func(matrix);
    
    matrix.CopyToHost();
    
    std::cout << "\nComplex CUDA Results:\n";
    std::cout << "  Row 0, Col 2: " << matrix[0][2] << "\n";
    std::cout << "  Row 0, Col 3: " << matrix[0][3] << "\n";
  }
  
  // ============================================================================
  // Generated CUDA Kernel (Example Output)
  // ============================================================================
  
  // For the first example, the system would generate approximately:
  /*
  
  __global__ void generated_kernel(
      double* matrix_0_data,  // matrixA
      size_t matrix_0_rows,
      size_t matrix_0_cols,
      double* matrix_1_data,  // matrixB
      size_t matrix_1_rows,
      size_t matrix_1_cols,
      size_t vector_length
  ) {
      int global_idx = blockIdx.x * blockDim.x + threadIdx.x;
      int group = global_idx / vector_length;
      int row_in_group = global_idx % vector_length;
      
      if (global_idx >= matrix_0_rows * vector_length) return;
      
      // Temporary variable
      double tmp_0;
      
      // First ForEachRow: tmp_0 = col(A,0) + col(B,2)
      tmp_0 = matrix_0_data[(group * matrix_0_cols + 0) * vector_length + row_in_group] +
              matrix_1_data[(group * matrix_1_cols + 2) * vector_length + row_in_group];
      
      // Second ForEachRow: col(A,1) = tmp_0
      matrix_0_data[(group * matrix_0_cols + 1) * vector_length + row_in_group] = tmp_0;
  }
  
  */
  
  return 0;
}

// ============================================================================
// Compilation:
//   nvcc -std=c++17 example_cuda_function.cpp -lnvrtc -lcuda
// ============================================================================
