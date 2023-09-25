#include <gtest/gtest.h>

#include <micm/process/jit_process_set.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>
#include <random>

#include "test_process_set_policy.hpp"

template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group2000VectorMatrix = micm::VectorMatrix<T, 2000>;
template<class T>
using Group3000VectorMatrix = micm::VectorMatrix<T, 3000>;
template<class T>
using Group4000VectorMatrix = micm::VectorMatrix<T, 4000>;

template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;

TEST(JitProcessSet, VectorMatrix)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  testProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix, micm::JitProcessSet<2>>(
      [&](const std::vector<micm::Process>& processes,
          const micm::State<Group2VectorMatrix>& state) -> micm::JitProcessSet<2> {
        return micm::JitProcessSet<2>{ jit.get(), processes, state };
      });
}

TEST(RandomJitProcessSet, VectorMatrix)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  testRandomSystem<Group2000VectorMatrix, micm::JitProcessSet<2000>>(
      2000,
      20,
      30,
      [&](const std::vector<micm::Process>& processes,
          const micm::State<Group2000VectorMatrix>& state) -> micm::JitProcessSet<2000> {
        return micm::JitProcessSet<2000>{ jit.get(), processes, state };
      });
  testRandomSystem<Group3000VectorMatrix, micm::JitProcessSet<3000>>(
      3000,
      50,
      40,
      [&](const std::vector<micm::Process>& processes,
          const micm::State<Group3000VectorMatrix>& state) -> micm::JitProcessSet<3000> {
        return micm::JitProcessSet<3000>{ jit.get(), processes, state };
      });
  testRandomSystem<Group3000VectorMatrix, micm::JitProcessSet<3000>>(
      3000,
      30,
      20,
      [&](const std::vector<micm::Process>& processes,
          const micm::State<Group3000VectorMatrix>& state) -> micm::JitProcessSet<3000> {
        return micm::JitProcessSet<3000>{ jit.get(), processes, state };
      });
  testRandomSystem<Group4000VectorMatrix, micm::JitProcessSet<4000>>(
      4000,
      100,
      80,
      [&](const std::vector<micm::Process>& processes,
          const micm::State<Group4000VectorMatrix>& state) -> micm::JitProcessSet<4000> {
        return micm::JitProcessSet<4000>{ jit.get(), processes, state };
      });
}