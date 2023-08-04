// Copyright (C) 2023 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <stdlib.h>

#include <iostream>
#include <memory>

#include "jit_compiler.hpp"
#include "llvm//IR/BasicBlock.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/LegacyPassManager.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Type.h"
#include "llvm/IR/Verifier.h"
#include "llvm/Support/TargetSelect.h"
#include "llvm/Target/TargetMachine.h"

namespace micm
{

  /// @brief Types used in JIT functions
  enum class JitType
  {
    Int,
    IntPtr,
    Float,
    FloatPtr,
    Double,
    DoublePtr,
    Bool,
    Void
  };

  /// @brief JIT function argument
  struct JitArgument
  {
    std::string name_;
    llvm::Type* type_;
    llvm::Value* arg_;
    llvm::AllocaInst* alloca_;
    llvm::Value* ptr_;
  };

  class JitFunctionBuilder;

  /// @brief A JIT-compiled function generator
  class JitFunction
  {
    std::string name_;
    std::shared_ptr<JitCompiler> compiler_;

   public:
    std::unique_ptr<llvm::LLVMContext> context_;
    std::unique_ptr<llvm::Module> module_;
    std::unique_ptr<llvm::IRBuilder<>> builder_;
    llvm::ExitOnError exit_on_error_;
    std::vector<JitArgument> arguments_;
    llvm::Function* function_;
    llvm::BasicBlock* entry_block_;

    JitFunction() = delete;

    friend class JitFunctionBuilder;
    static JitFunctionBuilder create(std::shared_ptr<JitCompiler> compiler);
    JitFunction(JitFunctionBuilder& function_builder);

    /// @brief Generates the function
    /// @return Resource tracker and function pointer
    ///
    /// This can only be called once.
    std::pair<llvm::orc::ResourceTrackerSP, llvm::JITTargetAddress> Generate();

   private:
    llvm::Type* GetType(JitType type);
    llvm::AllocaInst* CreateEntryBlockAlloca(llvm::Type* type, const std::string& var_name);
  };

  class JitFunctionBuilder
  {
    std::shared_ptr<JitCompiler> compiler_;
    std::string name_;
    std::vector<std::pair<std::string, JitType>> arguments_;
    JitType return_type_{ JitType::Void };
    friend class JitFunction;

   public:
    JitFunctionBuilder() = delete;
    JitFunctionBuilder(std::shared_ptr<JitCompiler> compiler);
    JitFunctionBuilder& name(std::string name);
    JitFunctionBuilder& arguments(const std::vector<std::pair<std::string, JitType>>& arguments);
    JitFunctionBuilder& return_type(JitType type);
  };

  inline JitFunctionBuilder JitFunction::create(std::shared_ptr<JitCompiler> compiler)
  {
    return JitFunctionBuilder{ compiler };
  }

  JitFunction::JitFunction(JitFunctionBuilder& function_builder)
      : name_(function_builder.name_),
        compiler_(function_builder.compiler_),
        context_(std::make_unique<llvm::LLVMContext>()),
        module_(std::make_unique<llvm::Module>(name_ + " module", *context_)),
        builder_(std::make_unique<llvm::IRBuilder<>>(*context_)),
        exit_on_error_(),
        arguments_(),
        function_(),
        entry_block_()
  {
    module_->setDataLayout(compiler_->GetDataLayout());

    // Prototype function
    std::vector<llvm::Type*> arg_types;
    for (const auto& pair : function_builder.arguments_)
    {
      llvm::Type* type = GetType(pair.second);
      arg_types.push_back(type);
      arguments_.push_back({ .name_ = pair.first, .type_ = type });
    }
    llvm::FunctionType* function_type = llvm::FunctionType::get(GetType(function_builder.return_type_), arg_types, false);
    function_ = llvm::Function::Create(function_type, llvm::Function::ExternalLinkage, name_, module_.get());
    llvm::Function::arg_iterator arg_iter = function_->arg_begin();
    for (auto& arg : arguments_)
    {
      arg.arg_ = arg_iter++;
      arg.arg_->setName(arg.name_);
    }

    // function body

    // set up entry block
    entry_block_ = llvm::BasicBlock::Create(*context_, "entry", function_);
    builder_->SetInsertPoint(entry_block_);

    // set up function argument variables
    for (auto& arg : arguments_)
      arg.alloca_ = CreateEntryBlockAlloca(arg.type_, arg.name_);
    for (auto& arg : arguments_)
      builder_->CreateStore(arg.arg_, arg.alloca_);
    for (auto& arg : arguments_)
      arg.ptr_ = builder_->CreateLoad(arg.type_, arg.alloca_);
  }

  std::pair<llvm::orc::ResourceTrackerSP, llvm::JITTargetAddress> JitFunction::Generate()
  {
    std::pair<llvm::orc::ResourceTrackerSP, llvm::JITTargetAddress> ret_val;
    verifyFunction(*function_);
    ret_val.first = compiler_->GetMainJITDylib().createResourceTracker();

    // Add the module to the JIT
    auto threadsafe_module = llvm::orc::ThreadSafeModule(std::move(module_), std::move(context_));
    exit_on_error_(compiler_->AddModule(std::move(threadsafe_module), ret_val.first));

    // Find the function
    auto expr_symbol = exit_on_error_(compiler_->Lookup(name_));
    ret_val.second = expr_symbol.getAddress();

    return ret_val;
  }

  inline llvm::Type* JitFunction::GetType(JitType type)
  {
    switch (type)
    {
      case JitType::Int: return llvm::Type::getInt64Ty(*context_);
      case JitType::IntPtr: return llvm::Type::getInt64Ty(*context_)->getPointerTo();
      case JitType::Float: return llvm::Type::getFloatTy(*context_);
      case JitType::FloatPtr: return llvm::Type::getFloatTy(*context_)->getPointerTo();
      case JitType::Double: return llvm::Type::getDoubleTy(*context_);
      case JitType::DoublePtr: return llvm::Type::getDoubleTy(*context_)->getPointerTo();
      case JitType::Bool: return llvm::Type::getInt1Ty(*context_);
      case JitType::Void: return llvm::Type::getVoidTy(*context_);
    }
    return llvm::Type::getVoidTy(*context_);
  }

  inline llvm::AllocaInst* JitFunction::CreateEntryBlockAlloca(llvm::Type* type, const std::string& var_name)
  {
    llvm::IRBuilder<> TmpB(&function_->getEntryBlock(), function_->getEntryBlock().begin());
    return TmpB.CreateAlloca(type, 0, var_name.c_str());
  }

  inline JitFunctionBuilder::JitFunctionBuilder(std::shared_ptr<JitCompiler> compiler)
      : compiler_(compiler){};

  inline JitFunctionBuilder& JitFunctionBuilder::name(std::string name)
  {
    name_ = name;
    return *this;
  }

  inline JitFunctionBuilder& JitFunctionBuilder::arguments(const std::vector<std::pair<std::string, JitType>>& arguments)
  {
    arguments_ = arguments;
    return *this;
  }

  inline JitFunctionBuilder& JitFunctionBuilder::return_type(JitType type)
  {
    return_type_ = type;
    return *this;
  }
}  // namespace micm