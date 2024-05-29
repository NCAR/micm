#include "analytical_policy.hpp"
#include "analytical_surface_rxn_policy.hpp"

#include <micm/solver/jit_rosenbrock.hpp>

#include <gtest/gtest.h>

template<class T>
using DefaultVectorMatrix = micm::VectorMatrix<T, 1>;

using DefaultSparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;

using DefaultJitRosenbrockSolver = micm::JitRosenbrockSolver<
    DefaultVectorMatrix,
    DefaultSparseVectorMatrix,
    micm::JitLinearSolver<1, DefaultSparseVectorMatrix, micm::JitLuDecomposition<1>>,
    micm::JitProcessSet<1>>;

TEST(AnalyticalExamplesJitRosenbrock, Troe)
{
  test_analytical_troe<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, TroeSuperStiffButAnalytical)
{
  test_analytical_stiff_troe<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, Photolysis)
{
  test_analytical_photolysis<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, PhotolysisSuperStiffButAnalytical)
{
  test_analytical_stiff_photolysis<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivation)
{
  test_analytical_ternary_chemical_activation<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, TernaryChemicalActivationSuperStiffButAnalytical)
{
  test_analytical_stiff_ternary_chemical_activation<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, Tunneling)
{
  test_analytical_tunneling<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, TunnelingSuperStiffButAnalytical)
{
  test_analytical_stiff_tunneling<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, Arrhenius)
{
  test_analytical_arrhenius<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, ArrheniusSuperStiffButAnalytical)
{
  test_analytical_stiff_arrhenius<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, Branched)
{
  test_analytical_branched<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, BranchedSuperStiffButAnalytical)
{
  test_analytical_stiff_branched<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, Robertson)
{
  test_analytical_robertson<DefaultJitRosenbrockSolver>();
}

TEST(AnalyticalExamplesJitRosenbrock, SurfaceRxn)
{
  test_analytical_surface_rxn<DefaultJitRosenbrockSolver>();
}
