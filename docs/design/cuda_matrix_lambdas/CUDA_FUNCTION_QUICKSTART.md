# Quick Start: Implementing CUDA Function() 

This guide provides concrete first steps to start implementing the CUDA JIT compilation system.

## Step 1: Set Up Development Environment

```bash
# Ensure CUDA toolkit is installed
nvcc --version

# Install NVRTC (usually comes with CUDA toolkit)
# Link libraries: -lnvrtc -lcuda

# Create feature branch
git checkout -b feature/cuda-function-jit
```

## Step 2: Create Minimal Test Case

Create `test/unit/cuda/util/test_cuda_function.cpp`:

```cpp
#include <gtest/gtest.h>
#include <micm/cuda/util/cuda_dense_matrix.hpp>

TEST(CudaFunction, MinimalAddition)
{
  // Start with the simplest possible case
  micm::CudaDenseMatrix<double, 1> matrix{3, 2, 0.0};
  
  // Initialize on host
  matrix[0][0] = 1.0;
  matrix[1][0] = 2.0;
  matrix[2][0] = 3.0;
  
  matrix.CopyToDevice();
  
  // Create a function that copies col(0) to col(1)
  // This avoids arithmetic operations initially
  auto func = micm::CudaDenseMatrix<double, 1>::Function(
      [](auto&& m) {
          m.ForEachRow(
              [&](auto a, auto& b) { b = a; },
              m.GetConstColumnView(0),
              m.GetColumnView(1));
      },
      matrix);
  
  func(matrix);
  
  matrix.CopyToHost();
  
  // Verify results
  EXPECT_EQ(matrix[0][1], 1.0);
  EXPECT_EQ(matrix[1][1], 2.0);
  EXPECT_EQ(matrix[2][1], 3.0);
}
```

## Step 3: Implement Core Expression Types

In `include/micm/cuda/util/cuda_function_expr.hpp`:

```cpp
namespace micm::cuda {

// Step 3a: Define expression operations
enum class ExprOp {
    LoadColumn,
    LoadTemp,
    Store,
    // Start with just these, add more later
};

// Step 3b: Define expression node
struct ExprNode {
    ExprOp op;
    int matrix_id = -1;
    int column_index = -1;
    int temp_var_id = -1;
    std::vector<std::shared_ptr<ExprNode>> children;
};

// Step 3c: Recording context to accumulate operations
struct RecordingContext {
    std::vector<std::shared_ptr<ExprNode>> operations;
    int num_matrices = 0;
    int num_temps = 0;
};

}  // namespace micm::cuda
```

## Step 4: Implement Recording Column View

In `include/micm/cuda/util/cuda_function_recording.hpp`:

```cpp
namespace micm::cuda {

class RecordingColumnView {
public:
    RecordingContext* ctx_;
    int matrix_id_;
    int column_index_;
    
    RecordingColumnView(RecordingContext* ctx, int mid, int cid)
        : ctx_(ctx), matrix_id_(mid), column_index_(cid) {}
    
    // Create an expression node representing this column
    std::shared_ptr<ExprNode> asExpr() const {
        auto node = std::make_shared<ExprNode>();
        node->op = ExprOp::LoadColumn;
        node->matrix_id = matrix_id_;
        node->column_index = column_index_;
        return node;
    }
};

}  // namespace micm::cuda
```

## Step 5: Implement Minimal Code Generator

In `include/micm/cuda/util/cuda_function_codegen.hpp`:

```cpp
namespace micm::cuda {

class CodeGenerator {
public:
    static std::string generateMinimalKernel(
        const RecordingContext& ctx,
        int src_col,
        int dst_col)
    {
        return R"(
__global__ void simple_copy_kernel(
    double* matrix_data,
    int num_rows,
    int num_cols,
    int src_col,
    int dst_col,
    int vector_length)
{
    int global_idx = blockIdx.x * blockDim.x + threadIdx.x;
    int group = global_idx / vector_length;
    int row_in_group = global_idx % vector_length;
    
    if (global_idx >= num_rows * vector_length) return;
    
    int src_index = (group * num_cols + src_col) * vector_length + row_in_group;
    int dst_index = (group * num_cols + dst_col) * vector_length + row_in_group;
    
    matrix_data[dst_index] = matrix_data[src_index];
}
)";
    }
};

}  // namespace micm::cuda
```

## Step 6: Implement NVRTC Compilation

In `include/micm/cuda/util/cuda_function_compiler.hpp`:

```cpp
#include <nvrtc.h>
#include <cuda.h>
#include <stdexcept>
#include <string>

namespace micm::cuda {

class NVRTCCompiler {
public:
    static CUfunction compile(const std::string& source) {
        // Create program
        nvrtcProgram prog;
        nvrtcCreateProgram(&prog, source.c_str(), "kernel.cu", 0, nullptr, nullptr);
        
        // Compile
        const char* opts[] = {"--gpu-architecture=compute_70"};
        nvrtcResult result = nvrtcCompileProgram(prog, 1, opts);
        
        if (result != NVRTC_SUCCESS) {
            size_t logSize;
            nvrtcGetProgramLogSize(prog, &logSize);
            std::string log(logSize, ' ');
            nvrtcGetProgramLog(prog, &log[0]);
            throw std::runtime_error("NVRTC compilation failed:\n" + log);
        }
        
        // Get PTX
        size_t ptxSize;
        nvrtcGetPTXSize(prog, &ptxSize);
        std::string ptx(ptxSize, ' ');
        nvrtcGetPTX(prog, &ptx[0]);
        
        nvrtcDestroyProgram(&prog);
        
        // Load module
        CUmodule module;
        cuModuleLoadDataEx(&module, ptx.c_str(), 0, nullptr, nullptr);
        
        // Get function
        CUfunction kernel;
        cuModuleGetFunction(&kernel, module, "simple_copy_kernel");
        
        return kernel;
    }
};

}  // namespace micm::cuda
```

## Step 7: Implement Function() Skeleton

In `include/micm/cuda/util/cuda_dense_matrix.hpp`, add:

```cpp
template<class T, std::size_t L>
class CudaDenseMatrix : public VectorMatrix<T, L> {
public:
    // ... existing methods ...
    
    template<typename Func, typename... Matrices>
    static auto Function(Func&& func, Matrices&... matrices) {
        // For MVP: Hard-code for single matrix, single copy operation
        
        // Step 1: Validate (simplified for now)
        static_assert(sizeof...(Matrices) == 1, "MVP supports single matrix only");
        
        // Step 2: Generate dummy kernel
        std::string kernel_source = cuda::CodeGenerator::generateMinimalKernel(
            /* ctx */ {}, 
            /* src_col */ 0, 
            /* dst_col */ 1);
        
        // Step 3: Compile
        CUfunction kernel = cuda::NVRTCCompiler::compile(kernel_source);
        
        // Step 4: Return callable
        return [kernel](Matrices&... invoke_matrices) {
            // Get first matrix
            auto& matrix = std::get<0>(std::tie(invoke_matrices...));
            
            // Prepare kernel parameters
            auto param = matrix.AsDeviceParam();
            void* args[] = {
                &param.d_data_,
                &matrix.NumRows(),
                &matrix.NumColumns(),
                /* src_col */ (void*)0,
                /* dst_col */ (void*)1,
                &L
            };
            
            // Launch
            int totalThreads = matrix.NumRows() * L;
            int blockSize = 256;
            int numBlocks = (totalThreads + blockSize - 1) / blockSize;
            
            cuLaunchKernel(
                kernel,
                numBlocks, 1, 1,
                blockSize, 1, 1,
                0, nullptr,
                args, nullptr);
            
            cuStreamSynchronize(nullptr);
        };
    }
};
```

## Step 8: Build and Test

```bash
# Add to CMakeLists.txt
find_package(CUDAToolkit REQUIRED)
target_link_libraries(your_target CUDA::nvrtc CUDA::cuda_driver)

# Build
mkdir build && cd build
cmake ..
make

# Run test
./test_cuda_function
```

## Step 9: Verify It Works

Expected output:
```
[==========] Running 1 test from 1 test suite.
[----------] 1 test from CudaFunction
[ RUN      ] CudaFunction.MinimalAddition
[       OK ] CudaFunction.MinimalAddition (250 ms)
[----------] 1 test from CudaFunction (250 ms total)
```

## Step 10: Iterate and Expand

Once the minimal test passes:

1. **Add arithmetic operations** (instead of just copy)
   ```cpp
   m.ForEachRow([&](auto a, auto& b) { b = a * 2.0; }, ...);
   ```

2. **Implement RecordingValue** with operator overloading
   ```cpp
   class RecordingValue {
       RecordingValue operator*(const RecordingValue& other);
   };
   ```

3. **Improve code generator** to handle expressions
   ```cpp
   std::string generateExpression(const ExprNode& node);
   ```

4. **Add more operations** (add, subtract, exp, log, etc.)

5. **Support multiple ForEachRow** calls

6. **Support multiple matrices**

7. **Add comprehensive tests**

## Common Issues and Solutions

### Issue 1: Linker errors with NVRTC
```bash
# Solution: Ensure proper linking
target_link_libraries(your_target PRIVATE CUDA::nvrtc CUDA::cuda_driver)
```

### Issue 2: Runtime compilation fails
```cpp
// Solution: Check CUDA toolkit installation
CUDA_ERROR_CHECK(cuInit(0));
```

### Issue 3: Kernel launch fails
```cpp
// Solution: Check for CUDA errors
CUresult result = cuLaunchKernel(...);
if (result != CUDA_SUCCESS) {
    const char* errStr;
    cuGetErrorString(result, &errStr);
    std::cerr << "Error: " << errStr << std::endl;
}
```

## Debugging Tips

1. **Print generated kernel:**
   ```cpp
   std::cout << "Generated kernel:\n" << kernel_source << std::endl;
   ```

2. **Use CUDA error checking:**
   ```cpp
   #define CUDA_CHECK(call) \
       do { \
           CUresult result = call; \
           if (result != CUDA_SUCCESS) { \
               const char* errStr; \
               cuGetErrorString(result, &errStr); \
               throw std::runtime_error(errStr); \
           } \
       } while(0)
   ```

3. **Test with simple kernels first:**
   Start with hardcoded kernels before automatic generation

4. **Validate on CPU first:**
   Compare every CUDA result with CPU implementation

## Success Criteria for MVP

- [ ] Test compiles without errors
- [ ] NVRTC compilation succeeds
- [ ] Kernel launches successfully
- [ ] Results match expected values
- [ ] No CUDA errors or memory leaks

Once MVP works, proceed to Phase 2 of the full implementation plan!
