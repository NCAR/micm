#include <gtest/gtest.h>

#include <micm/solver/jit_rosenbrock.hpp>

#include "analytical_policy.hpp"

template<class T>
using DefaultVectorMatrix = micm::VectorMatrix<T, 1>;

template<class T>
using DefaultSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;

using DefaultJitRosenbrockSolver = micm::JitRosenbrockSolver<
    DefaultVectorMatrix,
    DefaultSparseVectorMatrix,
    micm::JitLinearSolver<1, DefaultSparseVectorMatrix, micm::JitLuDecomposition<1>>>;

TEST(AnalyticalExamplesJitRosenbrock, Troe)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_troe<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}