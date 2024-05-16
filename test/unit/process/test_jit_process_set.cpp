#include "test_process_set_policy.hpp"

#include <micm/process/jit_process_set.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <random>

using Group2VectorMatrix = micm::VectorMatrix<double, 2>;
using Group200VectorMatrix = micm::VectorMatrix<double, 200>;
using Group300VectorMatrix = micm::VectorMatrix<double, 300>;
using Group400VectorMatrix = micm::VectorMatrix<double, 400>;

using Group2SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>;

TEST(JitProcessSet, VectorMatrix)
{
  auto jit{ micm::JitCompiler::Create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  testProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix, micm::JitProcessSet<2>>(
      [&](const std::vector<micm::Process>& processes,
          const micm::State<Group2VectorMatrix, Group2SparseVectorMatrix>& state) -> micm::JitProcessSet<2> {
        return micm::JitProcessSet<2>{ jit.get(), processes, state.variable_map_ };
      });
}

TEST(RandomJitProcessSet, VectorMatrix)
{
  auto jit{ micm::JitCompiler::Create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  testRandomSystem<Group200VectorMatrix, micm::StandardSparseMatrix, micm::JitProcessSet<200>>(
      200,
      20,
      30,
      [&](const std::vector<micm::Process>& processes,
          const micm::State<Group200VectorMatrix>& state) -> micm::JitProcessSet<200> {
        return micm::JitProcessSet<200>{ jit.get(), processes, state.variable_map_ };
      });
  testRandomSystem<Group300VectorMatrix, micm::StandardSparseMatrix, micm::JitProcessSet<300>>(
      300,
      50,
      40,
      [&](const std::vector<micm::Process>& processes,
          const micm::State<Group300VectorMatrix>& state) -> micm::JitProcessSet<300> {
        return micm::JitProcessSet<300>{ jit.get(), processes, state.variable_map_ };
      });
  testRandomSystem<Group300VectorMatrix, micm::StandardSparseMatrix, micm::JitProcessSet<300>>(
      300,
      30,
      20,
      [&](const std::vector<micm::Process>& processes,
          const micm::State<Group300VectorMatrix>& state) -> micm::JitProcessSet<300> {
        return micm::JitProcessSet<300>{ jit.get(), processes, state.variable_map_ };
      });
  testRandomSystem<Group400VectorMatrix, micm::StandardSparseMatrix, micm::JitProcessSet<400>>(
      400,
      100,
      80,
      [&](const std::vector<micm::Process>& processes,
          const micm::State<Group400VectorMatrix>& state) -> micm::JitProcessSet<400> {
        return micm::JitProcessSet<400>{ jit.get(), processes, state.variable_map_ };
      });
}