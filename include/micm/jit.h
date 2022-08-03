/* Copyright (C) 2022 National Center for Atmospheric Research,
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS),
 * and the U.S. Environmental Protection Agency (USEPA)
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef MICM_JIT_H
#define MICM_JIT_H

#include "llvm/ExecutionEngine/Orc/Core.h"
#include "llvm/ExecutionEngine/Orc/RTDyldObjectLinkingLayer.h"
#include "llvm/ExecutionEngine/Orc/IRCompileLayer.h"
#include "llvm/IR/DataLayout.h"
#include "llvm/ExecutionEngine/Orc/Mangling.h"
#include "llvm/ExecutionEngine/Orc/ThreadSafeModule.h"
#include "llvm/ExecutionEngine/Orc/JITTargetMachineBuilder.h"

class JIT {
    llvm::orc::ExecutionSession exec_session_;
    llvm::orc::RTDyldObjectLinkingLayer object_layer_;
    llvm::orc::IRCompileLayer compile_layer_;
    llvm::DataLayout::DataLayout data_layout_;
    llvm::orc::MangleAndInterner mangler_;
    llvm::orc::ThreadSafeContext context_;

  public:
    JIT(llvm::orc::JITTargetMachineBuilder builder, llvm::DataLayout::DataLayout data_layout)
};

#endif