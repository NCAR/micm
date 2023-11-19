#include <gtest/gtest.h>

#include <micm/solver/jit_rosenbrock.hpp>

#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

template<class T>
using DefaultVectorMatrix = micm::VectorMatrix<T, 1>;

template<class T>
using DefaultSparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;

using DefaultJitRosenbrockSolver = micm::JitRosenbrockSolver<
    DefaultVectorMatrix,
    DefaultSparseVectorMatrix,
    micm::JitLinearSolver<1, DefaultSparseVectorMatrix, micm::JitLuDecomposition<1>>,
    micm::JitProcessSet<1>>;

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

TEST(AnalyticalExamplesJitRosenbrock, TroeSuperStiffButAnalytical)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_stiff_troe<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, Photolysis)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_photolysis<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, PhotolysisSuperStiffButAnalytical)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_stiff_photolysis<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivation)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_ternary_chemical_activation<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_stiff_ternary_chemical_activation<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, Tunneling)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_tunneling<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, TunnelingSuperStiffButAnalytical)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_stiff_tunneling<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, Arrhenius)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_arrhenius<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, ArrheniusSuperStiffButAnalytical)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_stiff_arrhenius<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, Branched)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_branched<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, BranchedSuperStiffButAnalytical)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_stiff_branched<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, Robertson)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_robertson<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}

TEST(AnalyticalExamplesJitRosenbrock, SurfaceRxn)
{
  auto jit{ micm::JitCompiler::create() };
  if (auto err = jit.takeError())
  {
    llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
    EXPECT_TRUE(false);
  }
  test_analytical_surface_rxn<DefaultJitRosenbrockSolver>(
      [&](const micm::System& s, const std::vector<micm::Process>& p) -> DefaultJitRosenbrockSolver
      {
        return DefaultJitRosenbrockSolver{
          jit.get(), s, p, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
        };
      });
}
