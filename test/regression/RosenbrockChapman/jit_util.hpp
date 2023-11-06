#include <micm/solver/jit_rosenbrock.hpp>

#include "util.hpp"

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy,LinearSolverPolicy, ProcessSetPolicy> getTwoStageMultiCellJitChapmanSolver(
    const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  return micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy,LinearSolverPolicy, ProcessSetPolicy>(
      jit.get(),
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::move(processes),
      micm::RosenbrockSolverParameters::two_stage_rosenbrock_parameters(number_of_grid_cells));
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy,LinearSolverPolicy, ProcessSetPolicy> getThreeStageMultiCellJitChapmanSolver(
    const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  return micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy,LinearSolverPolicy, ProcessSetPolicy>(
      jit.get(),
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::move(processes),
      micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters(number_of_grid_cells));
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy,LinearSolverPolicy, ProcessSetPolicy> getFourStageMultiCellJitChapmanSolver(
    const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  return micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy,LinearSolverPolicy, ProcessSetPolicy>(
      jit.get(),
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::move(processes),
      micm::RosenbrockSolverParameters::four_stage_rosenbrock_parameters(number_of_grid_cells));
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy,LinearSolverPolicy, ProcessSetPolicy> getFourStageDAMultiCellJitChapmanSolver(
    const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  return micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy,LinearSolverPolicy, ProcessSetPolicy>(
      jit.get(),
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::move(processes),
      micm::RosenbrockSolverParameters::four_stage_differential_algebraic_rosenbrock_parameters(number_of_grid_cells));
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcessSetPolicy>
micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy,LinearSolverPolicy, ProcessSetPolicy> getSixStageDAMultiCellJitChapmanSolver(
    const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  return micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy,LinearSolverPolicy, ProcessSetPolicy>(
      jit.get(),
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::move(processes),
      micm::RosenbrockSolverParameters::six_stage_differential_algebraic_rosenbrock_parameters(number_of_grid_cells));
}