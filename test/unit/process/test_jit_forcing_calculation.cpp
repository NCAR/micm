// Copyright (C) 2023-2024 National Center for Atmospheric Research,
// SPDX-License-Identifier: Apache-2.0
//
// A comparison of the JITed forcing function with a hard-coded equivalent
// for a small chemical system.
//
// This code can be used as a test and a way to compare the IR of the
// JITed function with an optimized compilation of the hard-coded function

#include <gtest/gtest.h>

#include <chrono>
#include <micm/process/jit_process_set.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include "forcing_calculation.hpp"

#define NUM_ITERATIONS 100

template<class T>
using ForcingTestVectorMatrix = micm::VectorMatrix<T, NUM_GRID_CELLS>;
template<class T>
using ForcingTestSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<NUM_GRID_CELLS>>;

TEST(JitProcessSet, ForcingFunction)
{
  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");
  auto d = micm::Species("D");
  auto e = micm::Species("E");
  auto f = micm::Species("F");

  micm::Phase gas_phase{ std::vector<micm::Species>{ a, b, c, d, e, f } };

  micm::State<ForcingTestVectorMatrix, ForcingTestSparseVectorMatrix> state(
      micm::StateParameters{ .number_of_grid_cells_ = NUM_GRID_CELLS,
                             .number_of_rate_constants_ = 3,
                             .variable_names_{ "A", "B", "C", "D", "E", "F" } });

  micm::Process r1 = micm::Process::create().reactants({ a, b }).products({ micm::Yield{ d, 3.2 } }).phase(gas_phase);

  micm::Process r2 = micm::Process::create()
                         .reactants({ e, c })
                         .products({ micm::Yield{ a, 1.0 }, micm::Yield{ f, 2.0 } })
                         .phase(gas_phase);

  micm::Process r3 = micm::Process::create().reactants({ c, b }).products({ micm::Yield{ a, 1.0 } }).phase(gas_phase);

  std::vector<micm::Process> processes{ r1, r2, r3 };

  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  micm::JitProcessSet<NUM_GRID_CELLS> process_set{ jit.get(), processes, state.variable_map_ };

  ForcingTestVectorMatrix<double> rate_constants{ NUM_GRID_CELLS, 3 };
  ForcingTestVectorMatrix<double> forcing{ NUM_GRID_CELLS, 6, 0.0 };
  ForcingTestVectorMatrix<double> jit_forcing{ NUM_GRID_CELLS, 6, 0.0 };
  for (int i = 0; i < 3 * NUM_GRID_CELLS; ++i)
    rate_constants.AsVector()[i] = i * 2.0;
  for (int i = 0; i < 6 * NUM_GRID_CELLS; ++i)
    state.variables_.AsVector()[i] = i * 1.0;

  std::chrono::nanoseconds fixed_time{ 0 };
  std::chrono::nanoseconds jit_time{ 0 };
  for (std::size_t i = 0; i < NUM_ITERATIONS; ++i)
  {
    for (auto& elem : forcing.AsVector())
      elem = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
    calculate_forcing(rate_constants.AsVector().data(), state.variables_.AsVector().data(), forcing.AsVector().data());
    auto end = std::chrono::high_resolution_clock::now();
    fixed_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    for (auto& elem : jit_forcing.AsVector())
      elem = 0.0;
    start = std::chrono::high_resolution_clock::now();
    process_set.template AddForcingTerms<ForcingTestVectorMatrix>(rate_constants, state.variables_, jit_forcing);
    end = std::chrono::high_resolution_clock::now();
    jit_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    for (std::size_t i = 0; i < forcing.AsVector().size(); ++i)
    {
      EXPECT_EQ(forcing.AsVector()[i], jit_forcing.AsVector()[i]);
    }
  }
  std::cout << "Fixed forcing calculation time " << ((double)fixed_time.count()) / NUM_GRID_CELLS / NUM_ITERATIONS
            << " ns per grid cell" << std::endl;
  std::cout << "JITed forcing calculation time " << ((double)jit_time.count()) / NUM_GRID_CELLS / NUM_ITERATIONS
            << " ns per grid cell" << std::endl;
}