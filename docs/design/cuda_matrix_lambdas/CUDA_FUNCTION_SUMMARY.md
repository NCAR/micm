# CUDA Function() JIT Implementation - Summary and Next Steps

## What We're Building

A JIT compilation system for `CudaDenseMatrix::Function()` that:
1. Records matrix operations from user lambdas into an expression tree
2. Generates CUDA kernel source code from the expression tree
3. Compiles kernels using NVRTC at runtime
4. Caches compiled kernels for reuse
5. Launches kernels on the GPU with the same API as CPU matrices

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│ User Code                                                        │
│                                                                  │
│  auto func = CudaDenseMatrix::Function(                         │
│      [](auto&& m) {                                             │
│          auto tmp = m.GetRowVariable();                         │
│          m.ForEachRow([&](auto a, auto b, auto& t) {           │
│              t = a + b;                                         │
│          }, m.GetConstColumnView(0), ...);                      │
│      }, matrix);                                                │
│                                                                  │
│  func(matrix);  // Executes on GPU                             │
└─────────────────┬───────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────┐
│ Recording Phase                                                  │
│                                                                  │
│  - Execute lambda with RecordingGroupView proxies               │
│  - Operations return RecordingValue (not actual values)         │
│  - Build expression tree: Add(Load(col0), Load(col1))          │
│  - Track assignments, temporaries, column accesses              │
└─────────────────┬───────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────┐
│ Code Generation                                                  │
│                                                                  │
│  Walk expression tree → Generate CUDA kernel string:            │
│                                                                  │
│  __global__ void generated_kernel(...) {                        │
│      int idx = blockIdx.x * blockDim.x + threadIdx.x;          │
│      double tmp_0;                                              │
│      tmp_0 = matrix_data[...] + matrix_data[...];              │
│      matrix_data[...] = tmp_0;                                  │
│  }                                                              │
└─────────────────┬───────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────┐
│ NVRTC Compilation                                                │
│                                                                  │
│  - nvrtcCreateProgram(source)                                   │
│  - nvrtcCompileProgram(...)                                     │
│  - nvrtcGetPTX(...)                                             │
│  - cuModuleLoadDataEx(ptx)                                      │
│  - cuModuleGetFunction(...)                                     │
│  - Cache compiled kernel                                        │
└─────────────────┬───────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────┐
│ Kernel Execution                                                 │
│                                                                  │
│  - Set kernel parameters (matrix pointers, dimensions)          │
│  - cuLaunchKernel(...)                                          │
│  - cuStreamSynchronize(...)                                     │
└──────────────────────────────────────────────────────────────────┘
```

## Key Files Created

1. **[cuda_function_jit.md](cuda_function_jit.md)** - Complete design document
2. **[cuda_dense_matrix_function.hpp](../include/micm/cuda/util/cuda_dense_matrix_function.hpp)** - Expression templates and code generator
3. **[cuda_dense_matrix_function_impl.hpp](../include/micm/cuda/util/cuda_dense_matrix_function_impl.hpp)** - Function() implementation sketch
4. **[cuda_function_challenges.md](cuda_function_challenges.md)** - Implementation challenges and solutions
5. **[example_cuda_function.cpp](example_cuda_function.cpp)** - Usage examples

## Critical Implementation Details

### 1. Lambda Parameter Types

**Challenge**: C++ lambdas with `const double&` parameters can't directly accept `RecordingValue` objects.

**Solution**: Require `auto` parameters in nested lambdas:

```cpp
// CPU version (works today)
m.ForEachRow([&](const double& a, double& b) { b = a * 2; }, ...);

// CUDA version (proposed)
m.ForEachRow([&](auto a, auto& b) { b = a * 2; }, ...);
```

This allows the same lambda to work with both `double` (CPU) and `RecordingValue` (CUDA recording).

### 2. Expression Recording

Operations on `RecordingValue` objects build an expression tree:

```cpp
RecordingValue a = loadColumn(0);
RecordingValue b = loadColumn(1);
RecordingValue result = a + b;  // Creates Add(Load(0), Load(1))
temp = result;                   // Records assignment
```

### 3. Code Generation

Walk the expression tree to generate CUDA code:

```cpp
std::string generateExpression(ExprNode* node) {
    switch (node->op) {
        case ExprOp::LoadColumn:
            return "matrix_data[(group * cols + " + 
                   std::to_string(node->column_index) + 
                   ") * vector_length + row_in_group]";
        case ExprOp::Add:
            return "(" + generateExpression(node->children[0]) + 
                   " + " + generateExpression(node->children[1]) + ")";
        // ... other operations
    }
}
```

### 4. Kernel Caching

Hash the generated kernel source and cache the compiled function:

```cpp
std::unordered_map<std::string, CUfunction> kernel_cache_;

CUfunction getKernel(const std::string& source) {
    auto it = kernel_cache_.find(source);
    if (it != kernel_cache_.end()) {
        return it->second;  // Cache hit - no compilation
    }
    // Cache miss - compile and store
    CUfunction kernel = compileWithNVRTC(source);
    kernel_cache_[source] = kernel;
    return kernel;
}
```

## Implementation Phases

### Phase 1: MVP (Weeks 1-2)
**Goal**: Get one simple test case working end-to-end

- [ ] Implement basic expression nodes (Load, Add, Store)
- [ ] Implement RecordingValue with operator+
- [ ] Implement simple RecordingGroupView
- [ ] Implement basic code generator (single operation)
- [ ] Integrate NVRTC compilation
- [ ] Implement kernel cache
- [ ] Create test: `c = a + b` with single matrix

**Success Criteria**: One test passes comparing CUDA result to CPU

### Phase 2: Core Operations (Weeks 3-4)
**Goal**: Support all basic arithmetic and math functions

- [ ] Add operators: -, *, /
- [ ] Add math functions: exp, log, sqrt, pow, abs
- [ ] Support multiple temporary variables
- [ ] Support multiple ForEachRow calls (kernel fusion)
- [ ] Add comprehensive error checking

**Success Criteria**: testArrayFunction() passes with CUDA matrices

### Phase 3: Multi-Matrix (Week 5)
**Goal**: Support operations across multiple matrices

- [ ] Extend RecordingContext for multiple matrices
- [ ] Update code generator for multi-matrix kernels
- [ ] Add dimension validation
- [ ] Handle different matrix types (const vs non-const)

**Success Criteria**: testMultiMatrixArrayFunction() passes

### Phase 4: Testing & Optimization (Week 6)
**Goal**: Ensure correctness and performance

- [ ] Run all test_matrix_policy.hpp tests with CUDA
- [ ] Add CUDA-specific tests
- [ ] Profile kernel compilation time
- [ ] Profile kernel execution time
- [ ] Optimize generated code quality
- [ ] Add debug/logging infrastructure

**Success Criteria**: All tests pass, performance comparable to hand-written

### Phase 5: Polish (Week 7-8)
**Goal**: Production-ready code

- [ ] Comprehensive error messages
- [ ] Documentation and examples
- [ ] Handle edge cases (empty matrices, etc.)
- [ ] Support for different GPU architectures
- [ ] PTX caching to disk (optional)
- [ ] Integration with existing CUDA infrastructure

## Required Dependencies

1. **NVRTC** - NVIDIA Runtime Compilation library
   - Usually included with CUDA toolkit
   - Link with `-lnvrtc`

2. **CUDA Driver API** - For loading and launching kernels
   - Link with `-lcuda`

3. **C++17** - For template metaprogramming features
   - `if constexpr`, fold expressions, etc.

## API Considerations

### Option 1: Require 'auto' (Recommended)

```cpp
// Users must write:
m.ForEachRow([&](auto a, auto b, auto& c) { c = a + b; }, ...);
```

**Pros**: Simple, works for both CPU and CUDA
**Cons**: Less explicit type information

### Option 2: Separate API

```cpp
// Separate method for CUDA
m.CudaForEachRow([&](auto a, auto b, auto& c) { c = a + b; }, ...);
```

**Pros**: Clear separation, no confusion
**Cons**: Duplicate API, harder to switch between CPU/CUDA

### Option 3: Type Traits Detection

Detect at compile time whether using CUDA matrix and adapt accordingly.

**Pros**: Fully transparent
**Cons**: Complex template metaprogramming

**Recommendation**: Start with Option 1, consider Option 3 later if needed.

## Testing Strategy

### Unit Tests
```cpp
TEST(CudaFunction, BasicAddition) {
    CudaDenseMatrix<double, 1> m{3, 2, 0.0};
    // Initialize values
    auto func = CudaDenseMatrix<double, 1>::Function(...);
    func(m);
    m.CopyToHost();
    // Check results
}
```

### Comparison Tests
```cpp
template<template<class> class MatrixPolicy>
void testOperation() {
    // Same test for both VectorMatrix and CudaDenseMatrix
}

TEST(MatrixPolicy, CpuVersion) {
    testOperation<VectorMatrix>();
}

TEST(MatrixPolicy, CudaVersion) {
    testOperation<CudaDenseMatrix>();
}
```

### Performance Tests
```cpp
TEST(CudaFunction, CompilationTime) {
    auto start = std::chrono::high_resolution_clock::now();
    auto func = CudaDenseMatrix<double>::Function(...);
    auto compile_time = ... - start;
    EXPECT_LT(compile_time, std::chrono::milliseconds(1000));
}
```

## Risk Mitigation

### Risk 1: NVRTC Compilation Errors
**Mitigation**: Extensive testing, code generation validation, good error messages

### Risk 2: Performance Issues
**Mitigation**: Profile early, compare to hand-written kernels, optimize code generator

### Risk 3: Template Complexity
**Mitigation**: Incremental development, extensive comments, helper utilities

### Risk 4: Lambda Interception Limitations
**Mitigation**: Document requirements clearly, provide examples, consider alternatives

## Open Questions for Discussion

1. **GPU Architecture Support**: Target specific compute capabilities or runtime detection?
2. **Error Fallback**: If kernel compilation fails, fall back to CPU or throw?
3. **Debug Mode**: How verbose should debug output be by default?
4. **PTX Caching**: Cache to disk for faster startup on subsequent runs?
5. **Async Support**: Should kernel launches be async with CUDA streams?

## Next Steps

1. **Review this design** with the team
2. **Create feature branch** for implementation
3. **Start Phase 1 MVP** with minimal test case
4. **Weekly progress reviews** and iterative refinement
5. **Benchmark progress** against hand-written CUDA kernels

## Estimated Timeline

- **Phase 1 (MVP)**: 2 weeks
- **Phase 2 (Core)**: 2 weeks  
- **Phase 3 (Multi-matrix)**: 1 week
- **Phase 4 (Testing)**: 1 week
- **Phase 5 (Polish)**: 2 weeks

**Total**: ~8 weeks to production-ready implementation

## Success Metrics

1. ✅ All CPU tests pass with CUDA matrices
2. ✅ <500ms compilation time per unique kernel
3. ✅ Kernel execution time within 10% of hand-written CUDA
4. ✅ Less than 5% overhead from expression recording
5. ✅ Zero memory leaks or CUDA errors
6. ✅ Clear error messages for common mistakes

---

**This is an ambitious but achievable project.** The expression template + JIT approach gives you the flexibility you need while maintaining API compatibility with CPU matrices. The key is incremental development with frequent testing.
