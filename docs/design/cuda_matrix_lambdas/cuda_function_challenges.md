// CUDA Function() Implementation - Key Challenges and Solutions

## Challenge 1: Intercepting Lambda Parameter Bindings

**Problem**: When user writes:
```cpp
mA.ForEachRow([&](const double& a, const double& b, double& t) { t = a + b; },
              mA.GetConstColumnView(0), 
              mA.GetConstColumnView(1),
              tmp);
```

The lambda expects `double&` but we need to pass `RecordingValue` proxies.

**Solution Options**:

### Option A: Type Erasure with `auto`
Change the pattern to use `auto` parameters:
```cpp
mA.ForEachRow([&](auto a, auto b, auto& t) { t = a + b; },
              mA.GetConstColumnView(0), 
              mA.GetConstColumnView(1),
              tmp);
```
- ✅ Pros: Works with both recording and real execution
- ❌ Cons: Requires API change, less type safety

### Option B: Template Lambda Wrapper
Wrap the user's lambda with a template that adapts types:
```cpp
template<typename UserFunc>
auto adaptLambda(UserFunc&& func) {
    return [func](auto&&... args) {
        // Convert RecordingValue to double& when needed
        return func(convertArgs(args)...);
    };
}
```
- ✅ Pros: No API change needed
- ❌ Cons: Complex template metaprogramming

### Option C: Reflection/AST Parsing (C++20)
Use compile-time reflection to parse lambda:
```cpp
// Future C++ feature
constexpr auto ast = reflect(lambda);
// Analyze and recreate with different parameter types
```
- ✅ Pros: Clean, automatic
- ❌ Cons: Not yet standardized, limited compiler support

**Recommended: Option A** - Require `auto` parameters in lambdas used with CUDA matrices.

---

## Challenge 2: Capturing Assignment Operations

**Problem**: How do we know that `t = a + b` is an assignment to `t`?

In C++, `operator=` is special and can't be easily intercepted for recording.

**Solution**: Make RecordingValue track assignments

```cpp
class RecordingValue {
    std::shared_ptr<ExprNode> expr_;
    std::shared_ptr<ExprNode>* storage_ = nullptr;  // Where to store result
    
    RecordingValue& operator=(const RecordingValue& other) {
        if (storage_) {
            *storage_ = other.expr_;  // Record assignment
        }
        return *this;
    }
};
```

When creating a RecordingValue for an output parameter:
```cpp
RecordingValue createOutputProxy(RecordingRowVariable& var) {
    auto val = RecordingValue(ctx, var.asExpr());
    val.storage_ = &var.storage_location_;  // Track where assignments go
    return val;
}
```

---

## Challenge 3: Distinguishing Input vs Output Parameters

**Problem**: `ForEachRow` parameters can be:
- `const double&` → input (load value)
- `double&` → output (store result)

We need to know which is which to generate correct CUDA code.

**Solution**: Template parameter pack inspection

```cpp
template<typename Func, typename... Args>
void ForEachRow(Func&& func, Args&&... args) {
    // Use function_traits to extract lambda parameter types
    using FuncTraits = function_traits<std::decay_t<Func>>;
    
    // Check const-ness of each parameter type
    constexpr bool is_const[] = {
        std::is_const_v<std::remove_reference_t<
            typename FuncTraits::template arg<I>::type>>...
    };
    
    // Create appropriate proxies
    auto proxies = createProxies(is_const, std::forward<Args>(args)...);
    
    // Execute lambda
    std::apply(func, proxies);
}
```

---

## Challenge 4: Statement Ordering in Expression Tree

**Problem**: Operations like:
```cpp
tmp1 = a * b;
tmp2 = a + b;
c = tmp1 + tmp2;
```

Must be executed in order, but expression trees don't capture sequence.

**Solution**: Sequential statement list

Instead of one big expression tree, maintain a list of statements:

```cpp
struct ForEachRowCall {
    std::vector<Statement> statements;
    
    struct Statement {
        std::shared_ptr<ExprNode> target;  // What to assign to
        std::shared_ptr<ExprNode> value;   // What to assign
    };
};
```

Each assignment creates a new statement in order.

---

## Challenge 5: Multiple ForEachRow Calls

**Problem**: User can call ForEachRow multiple times:
```cpp
mA.ForEachRow([&](...) { tmp = a + b; }, ...);
mA.ForEachRow([&](...) { c = tmp * 2; }, ...);
```

Should this be one kernel or two?

**Solution**: Kernel fusion

Generate a single kernel with all operations:
```cuda
__global__ void kernel(...) {
    // First ForEachRow
    tmp = a + b;
    
    // Second ForEachRow  
    c = tmp * 2;
}
```

This maximizes performance by reducing kernel launch overhead.

---

## Challenge 6: Temporary Variable Allocation

**Problem**: Row variables need storage, but we're recording, not executing.

**Solution**: Virtual temporaries

During recording:
```cpp
RecordingRowVariable GetRowVariable() {
    int id = ctx_->allocate_temp();
    return RecordingRowVariable(ctx_, id);
}
```

In generated kernel:
```cuda
double tmp_0;  // Allocated for each temporary
double tmp_1;
```

---

## Challenge 7: Matrix Data Layout

**Problem**: VectorMatrix uses grouped layout:
```cpp
data_[(group * y_dim_ + column) * L + row_in_group]
```

CUDA kernel must use the same layout.

**Solution**: Template the indexing expression

```cpp
std::string generateColumnAccess(int matrix_id, int column) {
    return formatString(
        "matrix_%d_data[(group * matrix_%d_cols + %d) * vector_length + row_in_group]",
        matrix_id, matrix_id, column
    );
}
```

---

## Challenge 8: Error Handling and Debugging

**Problem**: Generated kernel might not compile, or might have bugs.

**Solution**: Multi-level diagnostics

1. **Pre-compilation validation**: Check expression tree for obvious errors
2. **NVRTC error reporting**: Capture and display compiler errors with context
3. **Kernel source output**: Allow printing generated code for inspection
4. **Runtime validation**: Check kernel launch errors

```cpp
#define MICM_CUDA_FUNCTION_DEBUG  // Enable debug output

auto func = CudaDenseMatrix<double>::Function([](auto&& m) {
    // Will print generated kernel code
    ...
});
```

---

## Challenge 9: Type Support (int vs double)

**Problem**: Matrices can be `CudaDenseMatrix<int>` or `<double>` or `<float>`.

**Solution**: Template the code generator

```cpp
template<typename T>
std::string CudaCodeGenerator::generateKernel(...) {
    std::string type_name = typeid(T).name();  // "double", "float", "int"
    
    code << "    " << type_name << " tmp_0;\n";
    // ...
}
```

---

## Challenge 10: Const Correctness

**Problem**: Some column views are const, others mutable.

**Solution**: Track mutability in RecordingColumnView

```cpp
class RecordingColumnView {
    bool is_mutable_;
    
    // During code generation:
    if (!is_mutable_) {
        // Can only read from this column
        // Generate load code
    } else {
        // Can write to this column
        // Generate store code
    }
};
```

---

## Implementation Roadmap

### Phase 1: Minimal Viable Product (MVP)
- [ ] Basic expression recording (Add, Sub, Mul, Div)
- [ ] Single matrix operations only
- [ ] Manual lambda with `auto` parameters
- [ ] Simple code generator
- [ ] NVRTC integration
- [ ] One test case working end-to-end

### Phase 2: Multi-Matrix Support
- [ ] Multiple input matrices
- [ ] Column views from different matrices
- [ ] Dimension validation

### Phase 3: Advanced Operations
- [ ] Math functions (exp, log, pow, sqrt)
- [ ] Multiple temporary variables
- [ ] Multiple ForEachRow calls (kernel fusion)

### Phase 4: Polish
- [ ] Error handling
- [ ] Debug output
- [ ] Performance optimization
- [ ] Documentation
- [ ] Comprehensive tests

---

## Testing Strategy

### Unit Tests
- Expression tree construction
- Code generation for each operation type
- NVRTC compilation
- Kernel cache hit/miss

### Integration Tests
- Compare CUDA results to CPU version
- All test cases from test_matrix_policy.hpp

### Performance Tests
- Compilation time measurement
- Kernel execution time vs hand-written
- Cache performance

---

## Open Questions

1. **Compiler Version**: What minimum CUDA version to target?
   - Recommend: CUDA 11.0+ (for better NVRTC features)

2. **GPU Architecture**: Support all architectures or specific ones?
   - Recommend: compute_70+ (Volta and newer)

3. **Error Recovery**: What happens if NVRTC compilation fails?
   - Fallback to CPU? Throw exception? Both with flag?

4. **Asynchronous Execution**: Support CUDA streams?
   - Phase 5 feature

5. **Multi-GPU**: How to handle multiple devices?
   - Phase 6 feature

---

## Alternative: Simpler Approach - String Injection

If expression templates prove too complex, consider allowing users to inject CUDA code directly:

```cpp
auto func = CudaDenseMatrix<double>::CudaFunction(R"(
    // CUDA code with placeholders
    tmp0 = COL(A, 0) + COL(B, 2);
    COL(A, 1) = tmp0;
)", matrixA, matrixB);
```

- ✅ Simpler implementation
- ❌ Less type-safe
- ❌ Different API from CPU version

**Not recommended** unless expression templates prove infeasible.
