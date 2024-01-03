// Copyright (C) 2023-2024 National Center for Atmospheric Research,
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <stdlib.h>

#include <cassert>
#include <iostream>
#include <memory>
#include <micm/util/random_string.hpp>

#include "jit_compiler.hpp"
#include "llvm/IR/Attributes.h"
#include "llvm/IR/BasicBlock.h"
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
    Int32,
    Int32Ptr,
    Int64,
    Int64Ptr,
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

  /// @brief JIT function loop
  struct JitLoop
  {
    std::string name_;
    llvm::BasicBlock* block_;
    llvm::PHINode* index_;
    llvm::Value* step_;
    llvm::Value* end_;
    llvm::BasicBlock* prior_block_;
    llvm::BasicBlock* after_block_;
  };

  class JitFunctionBuilder;

  /// @brief A JIT-compiled function generator
  ///
  /// An instance of this class can be used to build a single JIT function and includes
  /// some convenience functions for creating loops and operating on array elements
  class JitFunction
  {
    bool generated_ = false;
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

    /// @brief Get an LLVM type variable
    /// @param type Type to get
    /// @return LLVM type
    llvm::Type* GetType(JitType type);

    /// @brief Get a value from an array
    /// @param array_ptr Array pointer
    /// @param index Index in array to return value for
    /// @param type Data type of element
    /// @return Value of array element
    llvm::Value* GetArrayElement(JitArgument array_ptr, llvm::ArrayRef<llvm::Value*> index, JitType type);

    /// @brief Set the value in an array
    /// @param array_ptr Array pointer
    /// @param index Index in array to return value for
    /// @param type Data type of element
    /// @param value Value to set array element to
    void SetArrayElement(JitArgument array_ptr, llvm::ArrayRef<llvm::Value*> index, JitType type, llvm::Value* value);

    /// @brief Start a for loop
    /// @param name Label for the loop
    /// @param start Starting index
    /// @param end Ending index
    /// @param step Step size
    /// @return Loop reference
    JitLoop StartLoop(std::string name, int start, int end, int step);
    JitLoop StartLoop(std::string name, llvm::Value* start, llvm::Value* end, llvm::Value* step);

    /// @brief End a loop block
    /// @param loop Loop reference
    void EndLoop(JitLoop& loop);

   private:
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
      : generated_(false),
        name_(function_builder.name_ + generate_random_string()),
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
    for (unsigned int i = 0; i < arguments_.size(); ++i)
      function_->addParamAttr(i, llvm::Attribute::NoAlias);

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
    assert((!generated_) && "JIT Function already generated");
    std::pair<llvm::orc::ResourceTrackerSP, llvm::JITTargetAddress> ret_val;
    verifyFunction(*function_);
    ret_val.first = compiler_->GetMainJITDylib().createResourceTracker();

    // Add the module to the JIT
    auto threadsafe_module = llvm::orc::ThreadSafeModule(std::move(module_), std::move(context_));
    exit_on_error_(compiler_->AddModule(std::move(threadsafe_module), ret_val.first));

    // Find the function
    auto expr_symbol = exit_on_error_(compiler_->Lookup(name_));
    ret_val.second = expr_symbol.getAddress();

    generated_ = true;
    return ret_val;
  }

  inline llvm::Type* JitFunction::GetType(JitType type)
  {
    switch (type)
    {
      case JitType::Int32: return llvm::Type::getInt32Ty(*context_);
      case JitType::Int32Ptr: return llvm::Type::getInt32Ty(*context_)->getPointerTo();
      case JitType::Int64: return llvm::Type::getInt64Ty(*context_);
      case JitType::Int64Ptr: return llvm::Type::getInt64Ty(*context_)->getPointerTo();
      case JitType::Float: return llvm::Type::getFloatTy(*context_);
      case JitType::FloatPtr: return llvm::Type::getFloatTy(*context_)->getPointerTo();
      case JitType::Double: return llvm::Type::getDoubleTy(*context_);
      case JitType::DoublePtr: return llvm::Type::getDoubleTy(*context_)->getPointerTo();
      case JitType::Bool: return llvm::Type::getInt1Ty(*context_);
      case JitType::Void: return llvm::Type::getVoidTy(*context_);
    }
    return llvm::Type::getVoidTy(*context_);
  }

  llvm::Value* JitFunction::GetArrayElement(JitArgument array_ptr, llvm::ArrayRef<llvm::Value*> index, JitType type)
  {
    llvm::Value* elem = builder_->CreateGEP(GetType(type), array_ptr.ptr_, index, array_ptr.name_ + " get elem");
    return builder_->CreateLoad(GetType(type), elem, array_ptr.name_ + " load elem");
  }

  void
  JitFunction::SetArrayElement(JitArgument array_ptr, llvm::ArrayRef<llvm::Value*> index, JitType type, llvm::Value* value)
  {
    llvm::Value* elem = builder_->CreateGEP(GetType(type), array_ptr.ptr_, index, array_ptr.name_ + " set elem");
    builder_->CreateStore(value, elem);
  }

  JitLoop JitFunction::StartLoop(std::string name, int start, int end, int step = 1)
  {
    llvm::Value* start_val = llvm::ConstantInt::get(*context_, llvm::APInt(64, start));
    llvm::Value* step_val = llvm::ConstantInt::get(*context_, llvm::APInt(64, step));
    llvm::Value* end_val = llvm::ConstantInt::get(*context_, llvm::APInt(64, end));
    return StartLoop(name, start_val, end_val, step_val);
  }

  JitLoop JitFunction::StartLoop(std::string name, llvm::Value* start, llvm::Value* end, llvm::Value* step)
  {
    JitLoop loop;
    loop.name_ = name;
    loop.prior_block_ = builder_->GetInsertBlock();
    loop.block_ = llvm::BasicBlock::Create(*context_, name, function_);
    builder_->CreateBr(loop.block_);
    builder_->SetInsertPoint(loop.block_);
    loop.index_ = builder_->CreatePHI(GetType(JitType::Int64), 2, "i_" + name);
    loop.index_->addIncoming(start, loop.prior_block_);
    loop.step_ = step;
    loop.end_ = end;
    return loop;
  }

  void JitFunction::EndLoop(JitLoop& loop)
  {
    llvm::Value* nextIter = builder_->CreateNSWAdd(loop.index_, loop.step_, "next " + loop.name_);
    llvm::Value* atEnd = builder_->CreateICmpSGE(nextIter, loop.end_, "at end " + loop.name_);
    loop.after_block_ = llvm::BasicBlock::Create(*context_, "after " + loop.name_, function_);
    builder_->CreateCondBr(atEnd, loop.after_block_, loop.block_);
    builder_->SetInsertPoint(loop.after_block_);
    loop.index_->addIncoming(nextIter, loop.block_);
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