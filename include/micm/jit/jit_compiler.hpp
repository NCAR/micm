// Copyright (C) 2023-2025 University Corporation for Atmospheric Research-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Based on examples from the LLVM Project,
// under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#pragma once

#include <micm/util/error.hpp>

#include <llvm/ADT/StringRef.h>
#include <llvm/Analysis/LoopAccessAnalysis.h>
#include <llvm/ExecutionEngine/JITSymbol.h>
#include <llvm/ExecutionEngine/Orc/CompileUtils.h>
#include <llvm/ExecutionEngine/Orc/Core.h>
#include <llvm/ExecutionEngine/Orc/ExecutionUtils.h>
#include <llvm/ExecutionEngine/Orc/ExecutorProcessControl.h>
#include <llvm/ExecutionEngine/Orc/IRCompileLayer.h>
#include <llvm/ExecutionEngine/Orc/IRTransformLayer.h>
#include <llvm/ExecutionEngine/Orc/JITTargetMachineBuilder.h>
#include <llvm/ExecutionEngine/Orc/RTDyldObjectLinkingLayer.h>
#include <llvm/ExecutionEngine/SectionMemoryManager.h>
#include <llvm/IR/DataLayout.h>
#include <llvm/IR/LLVMContext.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/PassManager.h>
#include <llvm/Support/TargetSelect.h>
#include <llvm/Transforms/IPO.h>
#include <llvm/Transforms/InstCombine/InstCombine.h>
#include <llvm/Transforms/Scalar.h>
#include <llvm/Transforms/Scalar/GVN.h>
#include <llvm/Transforms/Utils.h>
#include <llvm/Transforms/Vectorize.h>

#include <memory>

enum class MicmJitErrc
{
  InvalidMatrix = MICM_JIT_ERROR_CODE_INVALID_MATRIX,
  MissingJitFunction = MICM_JIT_ERROR_CODE_MISSING_JIT_FUNCTION,
  FailedToBuild = MICM_JIT_ERROR_CODE_FAILED_TO_BUILD
};

namespace std
{
  template<>
  struct is_error_condition_enum<MicmJitErrc> : true_type
  {
  };
}  // namespace std

namespace
{
  class JitErrorCategory : public std::error_category
  {
   public:
    const char *name() const noexcept override
    {
      return MICM_ERROR_CATEGORY_JIT;
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmJitErrc>(ev))
      {
        case MicmJitErrc::InvalidMatrix:
          return "Invalid matrix for JIT compiled operation. Ensure matrix is Vector-ordered with vector dimension equal to "
                 "the nubmer of grid cells.";
        case MicmJitErrc::MissingJitFunction: return "Missing JIT-compiled function";
        default: return "Unknown error";
      }
    }
  };

  const JitErrorCategory JIT_ERROR{};
}  // namespace

inline std::error_code make_error_code(MicmJitErrc e)
{
  return { static_cast<int>(e), JIT_ERROR };
}

namespace micm
{

  // a singleton class
  class JitCompiler
  {
   private:
    std::unique_ptr<llvm::orc::ExecutionSession> execution_session_;

    llvm::DataLayout data_layout_;
    llvm::orc::MangleAndInterner mangle_;

    llvm::orc::RTDyldObjectLinkingLayer object_layer_;
    llvm::orc::IRCompileLayer compile_layer_;
    llvm::orc::IRTransformLayer optimize_layer_;

    llvm::orc::JITDylib &main_lib_;

   public:
    // Delete the copy constructor and assignment operator
    JitCompiler(const JitCompiler &) = delete;
    JitCompiler &operator=(const JitCompiler &) = delete;

    static JitCompiler &GetInstance()
    {
      static std::unique_ptr<JitCompiler> instance;
      if (!instance)
      {
        auto expectedInstance = Create();
        if (!expectedInstance)
        {
          throw std::system_error(make_error_code(MicmJitErrc::FailedToBuild));
        }
        instance = std::move(*expectedInstance);
      }
      return *instance;
    }

    ~JitCompiler()
    {
      if (auto Err = execution_session_->endSession())
        execution_session_->reportError(std::move(Err));
    }

    const llvm::DataLayout &GetDataLayout() const
    {
      return data_layout_;
    }

    llvm::orc::JITDylib &GetMainJITDylib()
    {
      return main_lib_;
    }

    llvm::Error AddModule(
        llvm::orc::ThreadSafeModule threadsafe_module,
        llvm::orc::ResourceTrackerSP resource_tracker = nullptr)
    {
      if (!resource_tracker)
        resource_tracker = main_lib_.getDefaultResourceTracker();
      return optimize_layer_.add(resource_tracker, std::move(threadsafe_module));
    }

    llvm::Expected<llvm::JITEvaluatedSymbol> Lookup(llvm::StringRef name)
    {
      return execution_session_->lookup({ &main_lib_ }, mangle_(name.str()));
    }

   private:
    static llvm::Expected<std::unique_ptr<JitCompiler>> Create()
    {
      llvm::InitializeNativeTarget();
      llvm::InitializeNativeTargetAsmPrinter();
      llvm::InitializeNativeTargetAsmParser();

      auto EPC = llvm::orc::SelfExecutorProcessControl::Create();
      if (!EPC)
        return EPC.takeError();

      auto execution_session = std::make_unique<llvm::orc::ExecutionSession>(std::move(*EPC));

      llvm::orc::JITTargetMachineBuilder machine_builder(execution_session->getExecutorProcessControl().getTargetTriple());

      auto data_layout = machine_builder.getDefaultDataLayoutForTarget();
      if (!data_layout)
        return data_layout.takeError();

      return llvm::Expected<std::unique_ptr<JitCompiler>>(std::unique_ptr<JitCompiler>(
          new JitCompiler(std::move(execution_session), std::move(machine_builder), std::move(*data_layout))));
    }

    JitCompiler(
        std::unique_ptr<llvm::orc::ExecutionSession> execution_session,
        llvm::orc::JITTargetMachineBuilder machine_builder,
        llvm::DataLayout data_layout)
        : execution_session_(std::move(execution_session)),
          data_layout_(std::move(data_layout)),
          mangle_(*this->execution_session_, this->data_layout_),
          object_layer_(*this->execution_session_, []() { return std::make_unique<llvm::SectionMemoryManager>(); }),
          compile_layer_(
              *this->execution_session_,
              object_layer_,
              std::make_unique<llvm::orc::ConcurrentIRCompiler>(std::move(machine_builder))),
          optimize_layer_(*this->execution_session_, compile_layer_, OptimizeModule),
          main_lib_(this->execution_session_->createBareJITDylib("<main>"))
    {
      main_lib_.addGenerator(
          llvm::cantFail(llvm::orc::DynamicLibrarySearchGenerator::GetForCurrentProcess(data_layout_.getGlobalPrefix())));
    }

    static llvm::Expected<llvm::orc::ThreadSafeModule> OptimizeModule(
        llvm::orc::ThreadSafeModule threadsafe_module,
        const llvm::orc::MaterializationResponsibility &responsibility)
    {
      threadsafe_module.withModuleDo(
          [](llvm::Module &module)
          {
            // Create a function pass manager.
            auto pass_manager = std::make_unique<llvm::legacy::FunctionPassManager>(&module);

            llvm::VectorizerParams::VectorizationFactor = 4;

            // Add some optimizations.
            // many from
            // https://www.intel.com/content/www/us/en/developer/articles/technical/optimize-llvm-code-data-analytics-vectorization.html#gs.6xkals
            // explanation of a few: https://seanforfun.github.io/llvm/2019/08/08/LLVMKaleidoscopeChap3.html
            // pass_manager->add(llvm::createFunctionInliningPass()); // causes segfault
            pass_manager->add(llvm::createLoopRotatePass());
            pass_manager->add(llvm::createLICMPass());
            pass_manager->add(llvm::createInstructionCombiningPass());
            pass_manager->add(llvm::createReassociatePass());
            pass_manager->add(llvm::createPromoteMemoryToRegisterPass());
            pass_manager->add(llvm::createGVNPass());
            pass_manager->add(llvm::createCFGSimplificationPass());
            pass_manager->add(llvm::createLoopVectorizePass());
            pass_manager->add(llvm::createSLPVectorizerPass());
            // pass_manager->add(llvm::createGlobalOptimizerPass()); // causes segfault
            pass_manager->doInitialization();

            // Run the optimizations over all functions in the module being added to
            // the JIT.
            for (auto &function : module)
            {
              pass_manager->run(function);
            }
          });

      return std::move(threadsafe_module);
    }
  };

}  // end namespace micm