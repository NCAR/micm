#include <gtest/gtest.h>

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/photolysis_rate_constant.hpp>
#include <micm/process/process.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/system.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <utility>
#include <vector>

using yields = std::pair<micm::Species, double>;

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy> getMultiCellChapmanSolver(const size_t number_of_grid_cells)
{
  auto o = micm::Species("O");
  auto o1d = micm::Species("O1D");
  auto o2 = micm::Species("O2");
  auto o3 = micm::Species("O3");
  auto m = micm::Species("M");
  auto ar = micm::Species("Ar");
  auto n2 = micm::Species("N2");
  auto h2o = micm::Species("H2O");
  auto co2 = micm::Species("CO2");

  micm::Phase gas_phase{ std::vector<micm::Species>{ m, ar, co2, h2o, n2, o1d, o, o2, o3 } };

  micm::Process r1 =
      micm::Process::create()
          .reactants({ o1d, n2 })
          .products({ yields(o, 1), yields(n2, 1) })
          .rate_constant(micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 2.15e-11, .C_ = 110 }))
          .phase(gas_phase);

  micm::Process r2 =
      micm::Process::create()
          .reactants({ o1d, o2 })
          .products({ yields(o, 1), yields(o2, 1) })
          .rate_constant(micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .C_ = 55 }))
          .phase(gas_phase);

  micm::Process r3 =
      micm::Process::create()
          .reactants({ o, o3 })
          .products({ yields(o2, 2) })
          .rate_constant(micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 8e-12, .C_ = -2060 }))
          .phase(gas_phase);

  micm::Process r4 =
      micm::Process::create()
          .reactants({ o, o2, m })
          .products({ yields(o3, 1), yields(m, 1) })
          .rate_constant(micm::ArrheniusRateConstant(micm::ArrheniusRateConstantParameters{ .A_ = 6.0e-34, .B_ = 2.4 }))
          .phase(gas_phase);

  micm::Process photo_1 = micm::Process::create()
                              .reactants({ o2 })
                              .products({ yields(o, 2) })
                              .rate_constant(micm::PhotolysisRateConstant())
                              .phase(gas_phase);

  micm::Process photo_2 = micm::Process::create()
                              .reactants({ o3 })
                              .products({ yields(o1d, 1), yields(o2, 1) })
                              .rate_constant(micm::PhotolysisRateConstant())
                              .phase(gas_phase);

  micm::Process photo_3 = micm::Process::create()
                              .reactants({ o3 })
                              .products({ yields(o, 1), yields(o2, 1) })
                              .rate_constant(micm::PhotolysisRateConstant())
                              .phase(gas_phase);

  return micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy>(
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::move(std::vector<micm::Process>{ photo_1, photo_2, photo_3, r1, r2, r3, r4 }),
      micm::three_stage_rosenbrock_parameters(number_of_grid_cells));
}
