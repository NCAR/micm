// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "../equilibrium_constraint_policy.hpp"

#include <micm/Kokkos.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

namespace
{
  auto MakeBuilder(const micm::RosenbrockSolverParameters& options)
  {
    return micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>(options);
  }
}  // namespace

TEST(EquilibriumIntegration, SetConstraintsAPIWorks)
{
  TestSetConstraintsAPIWorks(MakeBuilder);
}

TEST(EquilibriumIntegration, SetConstraintsAPIMultipleConstraints)
{
  TestSetConstraintsAPIMultipleConstraints(MakeBuilder);
}

TEST(EquilibriumIntegration, DAESolveWithConstraint)
{
  TestDAESolveWithConstraint(MakeBuilder);
}

TEST(EquilibriumIntegration, DAESolveWithConstraintAndReorderState)
{
  TestDAESolveWithConstraintAndReorderState(MakeBuilder);
}

TEST(EquilibriumIntegration, DAESolveWithTwoCoupledConstraints)
{
  TestDAESolveWithTwoCoupledConstraints(MakeBuilder);
}

TEST(EquilibriumIntegration, DAESolveWithNonUnitStoichiometry)
{
  TestDAESolveWithNonUnitStoichiometry(MakeBuilder);
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
