#include "util.hpp"

#include <micm/solver/jit_rosenbrock.hpp>

template<
    template<class>
    class MatrixPolicy,
    class SparseMatrixPolicy,
    class LinearSolverPolicy,
    class ProcessSetPolicy>
micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>
getTwoStageMultiCellJitChapmanSolver(const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  auto options = micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters(number_of_grid_cells);
  options.ignore_unused_species_ = true;

  return micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }), std::move(processes), options);
}

template<
    template<class>
    class MatrixPolicy,
    class SparseMatrixPolicy,
    class LinearSolverPolicy,
    class ProcessSetPolicy>
micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>
getThreeStageMultiCellJitChapmanSolver(const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  auto options = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters(number_of_grid_cells);
  options.ignore_unused_species_ = true;

  return micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }), std::move(processes), options);
}

template<
    template<class>
    class MatrixPolicy,
    class SparseMatrixPolicy,
    class LinearSolverPolicy,
    class ProcessSetPolicy>
micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>
getFourStageMultiCellJitChapmanSolver(const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  auto options = micm::RosenbrockSolverParameters::FourStageRosenbrockParameters(number_of_grid_cells);
  options.ignore_unused_species_ = true;

  return micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }), std::move(processes), options);
}

template<
    template<class>
    class MatrixPolicy,
    class SparseMatrixPolicy,
    class LinearSolverPolicy,
    class ProcessSetPolicy>
micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>
getFourStageDAMultiCellJitChapmanSolver(const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters(number_of_grid_cells);
  options.ignore_unused_species_ = true;

  return micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }), std::move(processes), options);
}

template<
    template<class>
    class MatrixPolicy,
    class SparseMatrixPolicy,
    class LinearSolverPolicy,
    class ProcessSetPolicy>
micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>
getSixStageDAMultiCellJitChapmanSolver(const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  auto options = micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters(number_of_grid_cells);
  options.ignore_unused_species_ = true;

  return micm::JitRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>( micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }), std::move(processes), options);
}