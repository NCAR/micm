#include "util.hpp"

#include <micm/jit/solver/jit_rosenbrock.hpp>
#include <micm/jit/solver/jit_solver_builder.hpp>
#include <micm/jit/solver/jit_solver_parameters.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

template<std::size_t L>
using JitBuilder = micm::JitSolverBuilder<micm::JitRosenbrockSolverParameters, L>;

template<std::size_t L>
auto getTwoStageMultiCellJitChapmanSolver(const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  return JitBuilder<L>(micm::RosenbrockSolverParameters::TwoStageRosenbrockParameters())
      .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
      .SetReactions(std::move(processes))
      .SetNumberOfGridCells(number_of_grid_cells)
      .SetIgnoreUnusedSpecies(true)
      .Build();
}

template<std::size_t L>
auto getThreeStageMultiCellJitChapmanSolver(const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  return JitBuilder<L>(micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
      .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
      .SetReactions(std::move(processes))
      .SetNumberOfGridCells(number_of_grid_cells)
      .SetIgnoreUnusedSpecies(true)
      .Build();
}

template<std::size_t L>
auto getFourStageMultiCellJitChapmanSolver(const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  return JitBuilder<L>(micm::RosenbrockSolverParameters::FourStageRosenbrockParameters())
      .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
      .SetReactions(std::move(processes))
      .SetNumberOfGridCells(number_of_grid_cells)
      .SetIgnoreUnusedSpecies(true)
      .Build();
}

template<std::size_t L>
auto getFourStageDAMultiCellJitChapmanSolver(const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  return JitBuilder<L>(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters())
      .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
      .SetReactions(std::move(processes))
      .SetNumberOfGridCells(number_of_grid_cells)
      .SetIgnoreUnusedSpecies(true)
      .Build();
}

template<std::size_t L>
auto getSixStageDAMultiCellJitChapmanSolver(const size_t number_of_grid_cells)
{
  micm::Phase gas_phase = createGasPhase();
  std::vector<micm::Process> processes = createProcesses(gas_phase);

  return JitBuilder<L>(micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters())
      .SetSystem(micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }))
      .SetReactions(std::move(processes))
      .SetNumberOfGridCells(number_of_grid_cells)
      .SetIgnoreUnusedSpecies(true)
      .Build();
}