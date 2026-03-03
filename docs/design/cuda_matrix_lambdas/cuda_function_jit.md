# CUDA Function() JIT Compilation Design

## Overview

This document describes the design for implementing `CudaDenseMatrix::Function()` with JIT compilation using NVRTC (NVIDIA Runtime Compilation). This enables users to write matrix operations using the same lambda syntax as CPU matrices, with automatic translation to compiled CUDA kernels.

## Goals

1. **API Compatibility**: Support the same lambda syntax as `VectorMatrix::Function()` and `Matrix::Function()`
2. **Performance**: Generate optimized CUDA kernels with minimal overhead
3. **Flexibility**: Support basic math operations on `ColumnView`s and `RowVariable`s
4. **Simplicity**: Start with essential operations, expand incrementally

## Architecture

### High-Level Flow

```
User Lambda → Expression Recording → Code Generation → NVRTC Compilation → Kernel Execution
```

### Components

#### 1. Expression Template System

Instead of performing actual operations, we create expression objects that represent the computation graph.

```cpp
// Expression types (not the actual values)
struct Expr {
    enum class Op { Add, Sub, Mul, Div, Exp, Log, Pow, Load, Store };
    Op operation;
    std::vector<Expr> children;
    // Metadata: variable ID, column index, constants, etc.
};
```

#### 2. Recording Proxy Objects

When `Function()` is called, matrix operations return **proxy objects** that record operations rather than executing them:

- `RecordingColumnView` - Records column access patterns
- `RecordingRowVariable` - Records temporary variable usage  
- `RecordingGroupView` - Orchestrates the recording

```cpp
// Example: User writes this
mA.ForEachRow([&](const double& a, const double& b, double& t) 
    { t = a + b; },
    mA.GetConstColumnView(0),
    mA.GetConstColumnView(1),
    tmp);

// We record: "tmp = col0 + col1" as an expression tree
```

#### 3. Code Generator

Walks the expression tree and generates CUDA kernel source code:

```cpp
std::string generateKernel(const std::vector<Operation>& ops) {
    return R"(
        __global__ void generated_kernel(
            double* matrix_data[],
            size_t num_rows,
            size_t num_cols[],
            size_t vector_length
        ) {
            int row = blockIdx.x * blockDim.x + threadIdx.x;
            if (row >= num_rows) return;
            
            // Generated operations
            double tmp0[MAX_VECTOR_SIZE];
            tmp0[tid] = matrix_data[0][row + col0 * num_rows] + 
                        matrix_data[1][row + col2 * num_rows];
            matrix_data[0][row + col1 * num_rows] = tmp0[tid];
        }
    )";
}
```

#### 4. NVRTC Compilation Pipeline

```cpp
class CudaKernelCache {
    // Compile kernel once, cache by expression signature
    CUfunction getOrCompile(const std::string& kernel_source);
    
    // Cache management
    std::unordered_map<std::string, CUfunction> cache_;
};
```

#### 5. Kernel Launcher

The returned callable object from `Function()` launches the compiled kernel:

```cpp
auto func = CudaDenseMatrix::Function([](auto&& m) { ... }, matrix);
func(matrix); // Launches CUDA kernel, not CPU loop
```

## Implementation Phases

### Phase 1: Core Infrastructure (MVP)
- [ ] Expression template base classes
- [ ] Recording proxy objects for ColumnView and RowVariable
- [ ] Simple operation recording (Add, Sub, Mul, Div)
- [ ] Basic code generator (single matrix, 1-2 operations)
- [ ] NVRTC integration and compilation
- [ ] Kernel cache
- [ ] Single matrix test case

### Phase 2: Extended Operations
- [ ] Math functions (exp, log, pow, sqrt, abs)
- [ ] Multiple temporary variables
- [ ] Multiple ForEachRow calls in sequence
- [ ] Const vs non-const column views
- [ ] Multi-matrix operations

### Phase 3: Optimization
- [ ] Kernel fusion (combine multiple ForEachRow into one kernel)
- [ ] Register usage optimization
- [ ] Shared memory for temporaries
- [ ] Automatic vectorization hints

### Phase 4: Robustness
- [ ] Error handling and diagnostics
- [ ] Compilation error reporting
- [ ] Runtime dimension validation
- [ ] PTX caching to disk
- [ ] Warmup/precompilation API

## Technical Details

### Expression Recording Strategy

Use **two-pass execution**:

1. **Recording Pass**: Execute lambda with proxy objects that build expression tree
   - No actual computation happens
   - Records: which columns accessed, what operations, what order
   
2. **Execution Pass**: Generate and run CUDA kernel
   - Translate expression tree to CUDA code
   - Compile with NVRTC
   - Launch kernel on device

### Handling Nested Lambdas

The challenge: `ForEachRow()` takes a lambda, but we can't introspect C++ lambdas.

**Solution**: Make `ForEachRow()` a template that captures the lambda's **type signature**:

```cpp
template<typename Func, typename... Args>
void RecordingGroupView::ForEachRow(Func&& func, Args&&... args) {
    // Extract lambda parameter types using template magic
    using Traits = function_traits<Func>;
    
    // Create "dummy" arguments that record operations
    auto result = func(createRecordingProxy(args)...);
    
    // The lambda body executed, recording operations
    // Now we have the expression tree
}
```

### Code Generation Pattern

For each `ForEachRow()` call, generate:

```cuda
// Allocate temporaries (if any)
double tmp_0, tmp_1;

// For each parameter, generate load/compute/store
if (is_const_column_view) {
    // Load: tmp_0 = matrix_data[col_idx]
}
if (is_row_variable) {
    // tmp_1 = <expression>
}
if (is_mutable_column_view) {
    // Store: matrix_data[col_idx] = <expression>
}
```

### NVRTC Compilation Example

```cpp
void compileKernel(const std::string& source) {
    nvrtcProgram prog;
    nvrtcCreateProgram(&prog, source.c_str(), "generated.cu", 0, nullptr, nullptr);
    
    const char* opts[] = {
        "--gpu-architecture=compute_70",
        "--std=c++17",
        "--use_fast_math"
    };
    
    nvrtcCompileProgram(prog, 3, opts);
    
    size_t ptxSize;
    nvrtcGetPTXSize(prog, &ptxSize);
    
    char* ptx = new char[ptxSize];
    nvrtcGetPTX(prog, ptx);
    
    // Load PTX and get function handle
    CUmodule module;
    cuModuleLoadDataEx(&module, ptx, 0, nullptr, nullptr);
    
    CUfunction kernel;
    cuModuleGetFunction(&kernel, module, "generated_kernel");
}
```

### Memory Layout Considerations

`VectorMatrix<T, L>` uses grouped layout:
- Data: `[(group * y_dim_ + column) * L + row_in_group]`
- Must generate correct indexing in CUDA kernel

### Performance Expectations

- **Compilation**: ~200-500ms per unique kernel (one-time cost)
- **Kernel Execution**: Near-optimal (equivalent to hand-written CUDA)
- **Cache Hit**: ~microseconds (lookup only)
- **Overall**: Amortized over many calls, negligible overhead

## API Example

```cpp
// User code - identical API to CPU version
CudaDenseMatrix<double, 4> matrixA{3, 2, 1.0};
CudaDenseMatrix<double, 4> matrixB{3, 3, 2.0};

// Copy initial data to device
matrixA.CopyToDevice();
matrixB.CopyToDevice();

// Create function (compiles kernel)
auto func = CudaDenseMatrix<double, 4>::Function(
    [](auto&& mA, auto&& mB) {
        auto tmp = mA.GetRowVariable();
        
        // First operation: tmp = mA.col(0) + mB.col(2)
        mA.ForEachRow([&](const double& a, const double& b, double& t) {
            t = a + b;
        }, mA.GetConstColumnView(0), mB.GetConstColumnView(2), tmp);
        
        // Second operation: mA.col(1) = tmp
        mA.ForEachRow([&](const double& t, double& c) {
            c = t;
        }, tmp, mA.GetColumnView(1));
    }, 
    matrixA, matrixB);

// Execute (launches kernel)
func(matrixA, matrixB);

// Copy results back
matrixA.CopyToHost();
```

## Alternative Approaches Considered

1. **Pre-defined Operations**: Too restrictive for chemistry simulations
2. **Extended Device Lambdas**: C++ lambda limitations with nested captures
3. **Macro-based DSL**: Poor user experience, non-idiomatic C++
4. **Expression Templates + JIT**: ✅ **Selected** - Best balance of flexibility and performance

## Dependencies

- **NVRTC**: NVIDIA Runtime Compilation library
- **CUDA Driver API**: For loading and launching compiled kernels
- **C++17**: Template metaprogramming features

## Testing Strategy

1. **Unit Tests**: Each operation type (add, mul, etc.)
2. **Integration Tests**: Multi-operation kernels
3. **Compatibility Tests**: Same results as CPU version
4. **Performance Tests**: Compare to hand-written CUDA kernels
5. **Compilation Tests**: Cache hit/miss behavior

## Open Questions

1. **Error Handling**: How to report NVRTC compilation errors to users?
2. **Debugging**: How to debug generated kernels?
3. **Portability**: Support for AMD (HIP) or other GPU backends?
4. **Conditional Logic**: Support for if/else in lambdas?
5. **Kernel Fusion**: Always fuse multiple ForEachRow, or make it optional?

## Future Extensions

- Support for reductions (sum, max, min across rows/columns)
- Custom user functions (via string injection)
- Automatic differentiation in generated kernels
- Multi-GPU support
- Async kernel launches with CUDA streams
- Graph API integration (CUDA 10+)

## References

- [NVRTC Documentation](https://docs.nvidia.com/cuda/nvrtc/)
- [Expression Templates in C++](https://en.wikipedia.org/wiki/Expression_templates)
- [Jitify Library](https://github.com/NVIDIA/jitify)
- Similar implementations: Eigen CUDA, Thrust, cuDF
