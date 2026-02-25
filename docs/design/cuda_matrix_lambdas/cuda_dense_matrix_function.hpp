// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/cuda/util/cuda_dense_matrix.hpp>
#include <cuda.h>
#include <nvrtc.h>

#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace micm
{
  namespace cuda
  {
    // ============================================================================
    // Expression Template System
    // ============================================================================

    /// @brief Types of operations that can be recorded
    enum class ExprOp
    {
      // Terminals
      LoadColumn,      // Load from matrix column
      LoadTemp,        // Load from temporary variable
      Constant,        // Literal constant value
      
      // Binary operations
      Add,
      Subtract,
      Multiply,
      Divide,
      
      // Unary operations
      Negate,
      
      // Math functions
      Exp,
      Log,
      Sqrt,
      Abs,
      Pow,             // Binary: base^exponent
      
      // Assignment
      Store            // Store to matrix column or temporary
    };

    /// @brief Represents a single expression node in the computation graph
    struct ExprNode
    {
      ExprOp op;
      
      // For terminals
      int source_matrix_id = -1;  // Which matrix (0, 1, 2, ...)
      int column_index = -1;       // Which column
      int temp_var_id = -1;        // Which temporary variable
      double constant_value = 0.0; // Constant value
      
      // For operations
      std::vector<std::shared_ptr<ExprNode>> children;
      
      ExprNode(ExprOp operation) : op(operation) {}
    };

    /// @brief Records a sequence of operations to be compiled into a kernel
    struct OperationSequence
    {
      struct ForEachRowCall
      {
        std::vector<std::shared_ptr<ExprNode>> statements;  // Statements in the lambda body
        std::vector<int> input_columns;                      // Column indices used as inputs
        std::vector<int> output_columns;                     // Column indices written to
        std::vector<int> temp_vars;                          // Temporary variable IDs used
      };
      
      std::vector<ForEachRowCall> foreach_calls;
      int num_temp_vars = 0;
      int num_matrices = 0;
      std::vector<std::pair<int, int>> matrix_dimensions;  // (rows, cols) for each matrix
    };

    // ============================================================================
    // Recording Proxy Objects
    // ============================================================================

    class RecordingContext;  // Forward declaration

    /// @brief Proxy for a column view that records accesses instead of reading/writing
    class RecordingColumnView
    {
     public:
      RecordingContext* ctx_;
      int matrix_id_;
      int column_index_;
      bool is_mutable_;
      
      RecordingColumnView(RecordingContext* ctx, int matrix_id, int col_idx, bool mutable_access)
          : ctx_(ctx), matrix_id_(matrix_id), column_index_(col_idx), is_mutable_(mutable_access)
      {
      }
      
      // Returns an expression node representing loading from this column
      std::shared_ptr<ExprNode> asExpr() const
      {
        auto node = std::make_shared<ExprNode>(ExprOp::LoadColumn);
        node->source_matrix_id = matrix_id_;
        node->column_index = column_index_;
        return node;
      }
    };

    /// @brief Proxy for a row variable (temporary storage)
    class RecordingRowVariable
    {
     public:
      RecordingContext* ctx_;
      int var_id_;
      
      RecordingRowVariable(RecordingContext* ctx, int id) : ctx_(ctx), var_id_(id) {}
      
      std::shared_ptr<ExprNode> asExpr() const
      {
        auto node = std::make_shared<ExprNode>(ExprOp::LoadTemp);
        node->temp_var_id = var_id_;
        return node;
      }
    };

    /// @brief Wrapper around expression nodes to enable operator overloading
    class RecordingValue
    {
     public:
      RecordingContext* ctx_;
      std::shared_ptr<ExprNode> expr_;
      
      RecordingValue(RecordingContext* ctx, std::shared_ptr<ExprNode> expr) 
          : ctx_(ctx), expr_(expr) {}
      
      // Arithmetic operators
      RecordingValue operator+(const RecordingValue& other) const;
      RecordingValue operator-(const RecordingValue& other) const;
      RecordingValue operator*(const RecordingValue& other) const;
      RecordingValue operator/(const RecordingValue& other) const;
      RecordingValue operator-() const;  // Unary negation
      
      // Math functions (implemented as friend functions below)
      friend RecordingValue exp(const RecordingValue& v);
      friend RecordingValue log(const RecordingValue& v);
      friend RecordingValue sqrt(const RecordingValue& v);
      friend RecordingValue abs(const RecordingValue& v);
      friend RecordingValue pow(const RecordingValue& base, const RecordingValue& exp);
    };

    /// @brief Context that accumulates recorded operations
    class RecordingContext
    {
     public:
      OperationSequence sequence_;
      std::vector<int> temp_var_stack_;  // Stack of available temporary variable IDs
      int next_temp_var_id_ = 0;
      
      // Allocate a new temporary variable
      RecordingRowVariable allocateTemp()
      {
        int id = next_temp_var_id_++;
        sequence_.num_temp_vars = std::max(sequence_.num_temp_vars, next_temp_var_id_);
        return RecordingRowVariable(this, id);
      }
      
      // Record a ForEachRow operation
      void recordForEachRow(const OperationSequence::ForEachRowCall& call)
      {
        sequence_.foreach_calls.push_back(call);
      }
      
      // Create binary operation expression
      std::shared_ptr<ExprNode> createBinaryOp(ExprOp op, std::shared_ptr<ExprNode> left, std::shared_ptr<ExprNode> right)
      {
        auto node = std::make_shared<ExprNode>(op);
        node->children = { left, right };
        return node;
      }
      
      // Create unary operation expression
      std::shared_ptr<ExprNode> createUnaryOp(ExprOp op, std::shared_ptr<ExprNode> operand)
      {
        auto node = std::make_shared<ExprNode>(op);
        node->children = { operand };
        return node;
      }
    };

    // Operator implementations
    inline RecordingValue RecordingValue::operator+(const RecordingValue& other) const
    {
      return RecordingValue(ctx_, ctx_->createBinaryOp(ExprOp::Add, expr_, other.expr_));
    }

    inline RecordingValue RecordingValue::operator-(const RecordingValue& other) const
    {
      return RecordingValue(ctx_, ctx_->createBinaryOp(ExprOp::Subtract, expr_, other.expr_));
    }

    inline RecordingValue RecordingValue::operator*(const RecordingValue& other) const
    {
      return RecordingValue(ctx_, ctx_->createBinaryOp(ExprOp::Multiply, expr_, other.expr_));
    }

    inline RecordingValue RecordingValue::operator/(const RecordingValue& other) const
    {
      return RecordingValue(ctx_, ctx_->createBinaryOp(ExprOp::Divide, expr_, other.expr_));
    }

    inline RecordingValue RecordingValue::operator-() const
    {
      return RecordingValue(ctx_, ctx_->createUnaryOp(ExprOp::Negate, expr_));
    }

    // Math functions
    inline RecordingValue exp(const RecordingValue& v)
    {
      return RecordingValue(v.ctx_, v.ctx_->createUnaryOp(ExprOp::Exp, v.expr_));
    }

    inline RecordingValue log(const RecordingValue& v)
    {
      return RecordingValue(v.ctx_, v.ctx_->createUnaryOp(ExprOp::Log, v.expr_));
    }

    inline RecordingValue sqrt(const RecordingValue& v)
    {
      return RecordingValue(v.ctx_, v.ctx_->createUnaryOp(ExprOp::Sqrt, v.expr_));
    }

    inline RecordingValue abs(const RecordingValue& v)
    {
      return RecordingValue(v.ctx_, v.ctx_->createUnaryOp(ExprOp::Abs, v.expr_));
    }

    inline RecordingValue pow(const RecordingValue& base, const RecordingValue& exponent)
    {
      return RecordingValue(base.ctx_, base.ctx_->createBinaryOp(ExprOp::Pow, base.expr_, exponent.expr_));
    }

    // ============================================================================
    // Recording Group View - Mimics normal GroupView but records operations
    // ============================================================================

    template<std::size_t L>
    class RecordingGroupView
    {
     public:
      RecordingContext* ctx_;
      int matrix_id_;
      
      RecordingGroupView(RecordingContext* ctx, int matrix_id) 
          : ctx_(ctx), matrix_id_(matrix_id) {}
      
      // Get a const column view (for reading)
      RecordingColumnView GetConstColumnView(std::size_t column_index) const
      {
        return RecordingColumnView(ctx_, matrix_id_, column_index, false);
      }
      
      // Get a mutable column view (for writing)
      RecordingColumnView GetColumnView(std::size_t column_index)
      {
        return RecordingColumnView(ctx_, matrix_id_, column_index, true);
      }
      
      // Get a temporary row variable
      RecordingRowVariable GetRowVariable()
      {
        return ctx_->allocateTemp();
      }
      
      // ForEachRow - This is where the magic happens
      // We execute the lambda with RecordingValue proxies to capture the expression
      template<typename Func, typename... Args>
      void ForEachRow(Func&& func, Args&&... args)
      {
        // TODO: This is the complex part - we need to:
        // 1. Convert args (ColumnViews, RowVariables) to RecordingValue proxies
        // 2. Call func with these proxies
        // 3. Record the operations that happen inside func
        // 4. Extract assignment targets
        
        // For now, placeholder - full implementation requires careful template metaprogramming
        OperationSequence::ForEachRowCall call;
        // ... recording logic ...
        ctx_->recordForEachRow(call);
      }
      
      std::size_t NumRows() const { return 0; }  // Placeholder
      std::size_t NumColumns() const { return 0; }  // Placeholder
    };

    // ============================================================================
    // CUDA Code Generator
    // ============================================================================

    class CudaCodeGenerator
    {
     public:
      /// @brief Generate CUDA kernel source code from operation sequence
      static std::string generateKernel(const OperationSequence& ops, const std::string& kernel_name = "generated_kernel")
      {
        std::ostringstream code;
        
        // Kernel signature
        code << "__global__ void " << kernel_name << "(\n";
        for (int i = 0; i < ops.num_matrices; ++i)
        {
          code << "    double* matrix_" << i << "_data,\n";
          code << "    size_t matrix_" << i << "_rows,\n";
          code << "    size_t matrix_" << i << "_cols,\n";
        }
        code << "    size_t vector_length\n";
        code << ") {\n";
        
        // Thread indexing
        code << "    int global_idx = blockIdx.x * blockDim.x + threadIdx.x;\n";
        code << "    int group = global_idx / vector_length;\n";
        code << "    int row_in_group = global_idx % vector_length;\n";
        code << "    \n";
        code << "    if (global_idx >= matrix_0_rows * vector_length) return;\n";
        code << "    \n";
        
        // Declare temporary variables
        if (ops.num_temp_vars > 0)
        {
          for (int i = 0; i < ops.num_temp_vars; ++i)
          {
            code << "    double tmp_" << i << ";\n";
          }
          code << "    \n";
        }
        
        // Generate code for each ForEachRow call
        for (const auto& call : ops.foreach_calls)
        {
          code << "    // ForEachRow operation\n";
          for (const auto& stmt : call.statements)
          {
            code << "    " << generateExpression(stmt) << ";\n";
          }
          code << "    \n";
        }
        
        code << "}\n";
        
        return code.str();
      }
      
     private:
      /// @brief Generate C code for an expression node
      static std::string generateExpression(const std::shared_ptr<ExprNode>& node)
      {
        switch (node->op)
        {
          case ExprOp::LoadColumn:
            // VectorMatrix layout: data_[(group * y_dim_ + column) * L + row_in_group]
            return "matrix_" + std::to_string(node->source_matrix_id) + "_data[" +
                   "(group * matrix_" + std::to_string(node->source_matrix_id) + "_cols + " +
                   std::to_string(node->column_index) + ") * vector_length + row_in_group]";
          
          case ExprOp::LoadTemp:
            return "tmp_" + std::to_string(node->temp_var_id);
          
          case ExprOp::Constant:
            return std::to_string(node->constant_value);
          
          case ExprOp::Add:
            return "(" + generateExpression(node->children[0]) + " + " + generateExpression(node->children[1]) + ")";
          
          case ExprOp::Subtract:
            return "(" + generateExpression(node->children[0]) + " - " + generateExpression(node->children[1]) + ")";
          
          case ExprOp::Multiply:
            return "(" + generateExpression(node->children[0]) + " * " + generateExpression(node->children[1]) + ")";
          
          case ExprOp::Divide:
            return "(" + generateExpression(node->children[0]) + " / " + generateExpression(node->children[1]) + ")";
          
          case ExprOp::Negate:
            return "(-" + generateExpression(node->children[0]) + ")";
          
          case ExprOp::Exp:
            return "exp(" + generateExpression(node->children[0]) + ")";
          
          case ExprOp::Log:
            return "log(" + generateExpression(node->children[0]) + ")";
          
          case ExprOp::Sqrt:
            return "sqrt(" + generateExpression(node->children[0]) + ")";
          
          case ExprOp::Abs:
            return "abs(" + generateExpression(node->children[0]) + ")";
          
          case ExprOp::Pow:
            return "pow(" + generateExpression(node->children[0]) + ", " + generateExpression(node->children[1]) + ")";
          
          default:
            return "/* unknown operation */";
        }
      }
    };

    // ============================================================================
    // NVRTC Kernel Compiler and Cache
    // ============================================================================

    class KernelCompiler
    {
     public:
      /// @brief Compile CUDA kernel from source code using NVRTC
      static CUfunction compile(const std::string& source, const std::string& kernel_name)
      {
        nvrtcProgram prog;
        nvrtcResult result = nvrtcCreateProgram(
            &prog,
            source.c_str(),
            "generated_kernel.cu",
            0,
            nullptr,
            nullptr);
        
        if (result != NVRTC_SUCCESS)
        {
          throw std::runtime_error("Failed to create NVRTC program");
        }
        
        // Compilation options
        const char* opts[] = {
            "--gpu-architecture=compute_70",  // Adjust based on target GPU
            "--std=c++17",
            "--use_fast_math",
            "--fmad=true"
        };
        
        result = nvrtcCompileProgram(prog, 4, opts);
        
        if (result != NVRTC_SUCCESS)
        {
          // Get compilation log
          size_t logSize;
          nvrtcGetProgramLogSize(prog, &logSize);
          std::vector<char> log(logSize);
          nvrtcGetProgramLog(prog, log.data());
          
          nvrtcDestroyProgram(&prog);
          throw std::runtime_error(std::string("NVRTC compilation failed:\n") + log.data());
        }
        
        // Get PTX
        size_t ptxSize;
        nvrtcGetPTXSize(prog, &ptxSize);
        std::vector<char> ptx(ptxSize);
        nvrtcGetPTX(prog, ptx.data());
        
        nvrtcDestroyProgram(&prog);
        
        // Load module and get function
        CUmodule module;
        CUresult cuResult = cuModuleLoadDataEx(&module, ptx.data(), 0, nullptr, nullptr);
        if (cuResult != CUDA_SUCCESS)
        {
          throw std::runtime_error("Failed to load CUDA module");
        }
        
        CUfunction kernel;
        cuResult = cuModuleGetFunction(&kernel, module, kernel_name.c_str());
        if (cuResult != CUDA_SUCCESS)
        {
          throw std::runtime_error("Failed to get kernel function");
        }
        
        return kernel;
      }
    };

    /// @brief Cache compiled kernels to avoid recompilation
    class KernelCache
    {
     private:
      std::unordered_map<std::string, CUfunction> cache_;
      
     public:
      /// @brief Get compiled kernel, compiling if not in cache
      CUfunction getKernel(const std::string& source, const std::string& kernel_name = "generated_kernel")
      {
        // Use source code as cache key (could use hash for efficiency)
        auto it = cache_.find(source);
        if (it != cache_.end())
        {
          return it->second;
        }
        
        // Not in cache, compile it
        CUfunction kernel = KernelCompiler::compile(source, kernel_name);
        cache_[source] = kernel;
        return kernel;
      }
      
      static KernelCache& instance()
      {
        static KernelCache cache;
        return cache;
      }
    };

  }  // namespace cuda
}  // namespace micm
