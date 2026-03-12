// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/linear_constraint.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <system_error>
#include <vector>

using namespace micm;

TEST(LinearConstraint, SimpleConservation)
{
  // Test: A + B = 1.0 (total concentration conservation)
  // Constraint: G = [A] + [B] - 1.0 = 0

  LinearConstraint constraint(
      "A_B_conservation",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      1.0);

  EXPECT_EQ(constraint.name_, "A_B_conservation");
  EXPECT_EQ(constraint.constant_, 1.0);
  EXPECT_EQ(constraint.species_dependencies_.size(), 2);
  EXPECT_EQ(constraint.species_dependencies_[0], "A");
  EXPECT_EQ(constraint.species_dependencies_[1], "B");
  EXPECT_EQ(constraint.terms_.size(), 2);

  // Test when constraint is satisfied: [A] = 0.3, [B] = 0.7
  // Residual = 0.3 + 0.7 - 1.0 = 0
  std::vector<double> concentrations = { 0.3, 0.7 };
  std::vector<std::size_t> indices = { 0, 1 };

  double residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 0.0, 1e-10);

  // Test when constraint is not satisfied: [A] = 0.5, [B] = 0.6
  // Residual = 0.5 + 0.6 - 1.0 = 0.1
  concentrations = { 0.5, 0.6 };
  residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 0.1, 1e-10);
}

TEST(LinearConstraint, Jacobian)
{
  // Test Jacobian for A + B = 1.0
  // G = [A] + [B] - 1.0
  // dG/d[A] = 1.0
  // dG/d[B] = 1.0

  LinearConstraint constraint(
      "A_B_conservation",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      1.0);

  std::vector<double> concentrations = { 0.3, 0.7 };
  std::vector<std::size_t> indices = { 0, 1 };
  std::vector<double> jacobian(2);

  constraint.Jacobian(concentrations.data(), indices.data(), jacobian.data());

  EXPECT_NEAR(jacobian[0], 1.0, 1e-10);  // dG/d[A] = 1.0
  EXPECT_NEAR(jacobian[1], 1.0, 1e-10);  // dG/d[B] = 1.0
}

TEST(LinearConstraint, WeightedSum)
{
  // Test: 2*A + 3*B - C = 5.0
  // Constraint: G = 2*[A] + 3*[B] - [C] - 5.0 = 0

  LinearConstraint constraint(
      "weighted_sum",
      std::vector<StoichSpecies>{
          StoichSpecies(Species("A"), 2.0),
          StoichSpecies(Species("B"), 3.0),
          StoichSpecies(Species("C"), -1.0) },
      5.0);

  EXPECT_EQ(constraint.species_dependencies_.size(), 3);
  EXPECT_EQ(constraint.constant_, 5.0);

  // Test when constraint is satisfied: [A] = 1.0, [B] = 2.0, [C] = 3.0
  // Residual = 2*1.0 + 3*2.0 - 3.0 - 5.0 = 2 + 6 - 3 - 5 = 0
  std::vector<double> concentrations = { 1.0, 2.0, 3.0 };
  std::vector<std::size_t> indices = { 0, 1, 2 };

  double residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 0.0, 1e-10);

  // Test Jacobian
  std::vector<double> jacobian(3);
  constraint.Jacobian(concentrations.data(), indices.data(), jacobian.data());

  EXPECT_NEAR(jacobian[0], 2.0, 1e-10);   // dG/d[A] = 2.0
  EXPECT_NEAR(jacobian[1], 3.0, 1e-10);   // dG/d[B] = 3.0
  EXPECT_NEAR(jacobian[2], -1.0, 1e-10);  // dG/d[C] = -1.0
}

TEST(LinearConstraint, ThreeSpeciesConservation)
{
  // Test: A + B + C = 10.0
  // Constraint: G = [A] + [B] + [C] - 10.0 = 0

  LinearConstraint constraint(
      "ABC_total",
      std::vector<StoichSpecies>{
          StoichSpecies(Species("A"), 1.0),
          StoichSpecies(Species("B"), 1.0),
          StoichSpecies(Species("C"), 1.0) },
      10.0);

  EXPECT_EQ(constraint.species_dependencies_.size(), 3);
  EXPECT_EQ(constraint.constant_, 10.0);

  // Test when constraint is satisfied: [A] = 2.0, [B] = 3.0, [C] = 5.0
  std::vector<double> concentrations = { 2.0, 3.0, 5.0 };
  std::vector<std::size_t> indices = { 0, 1, 2 };

  double residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 0.0, 1e-10);

  // Test when constraint is violated: [A] = 2.0, [B] = 3.0, [C] = 6.0
  // Residual = 2.0 + 3.0 + 6.0 - 10.0 = 1.0
  concentrations = { 2.0, 3.0, 6.0 };
  residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 1.0, 1e-10);
}

TEST(LinearConstraint, AlgebraicSpecies)
{
  // Test that AlgebraicSpecies returns the last species in terms list
  LinearConstraint constraint(
      "A_B_C_conservation",
      std::vector<StoichSpecies>{
          StoichSpecies(Species("A"), 1.0),
          StoichSpecies(Species("B"), 1.0),
          StoichSpecies(Species("C"), 1.0) },
      10.0);

  EXPECT_EQ(constraint.AlgebraicSpecies(), "C");
}

TEST(LinearConstraint, ZeroConstant)
{
  // Test: A - B = 0 (species balance)
  // Constraint: G = [A] - [B] = 0

  LinearConstraint constraint(
      "A_equals_B",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), -1.0) },
      0.0);

  EXPECT_EQ(constraint.constant_, 0.0);

  // Test when constraint is satisfied: [A] = 0.5, [B] = 0.5
  std::vector<double> concentrations = { 0.5, 0.5 };
  std::vector<std::size_t> indices = { 0, 1 };

  double residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 0.0, 1e-10);

  // Test when constraint is violated: [A] = 0.6, [B] = 0.5
  // Residual = 0.6 - 0.5 = 0.1
  concentrations = { 0.6, 0.5 };
  residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 0.1, 1e-10);
}

TEST(LinearConstraint, FractionalCoefficients)
{
  // Test: 0.5*A + 1.5*B = 2.0
  // Constraint: G = 0.5*[A] + 1.5*[B] - 2.0 = 0

  LinearConstraint constraint(
      "fractional_conservation",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 0.5), StoichSpecies(Species("B"), 1.5) },
      2.0);

  // Test when constraint is satisfied: [A] = 2.0, [B] = 1.0
  // 0.5*2.0 + 1.5*1.0 = 1.0 + 1.5 = 2.5... wait that's 2.5
  // Let's use [A] = 1.0, [B] = 1.0
  // 0.5*1.0 + 1.5*1.0 = 0.5 + 1.5 = 2.0 ✓
  std::vector<double> concentrations = { 1.0, 1.0 };
  std::vector<std::size_t> indices = { 0, 1 };

  double residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 0.0, 1e-10);

  // Test Jacobian
  std::vector<double> jacobian(2);
  constraint.Jacobian(concentrations.data(), indices.data(), jacobian.data());

  EXPECT_NEAR(jacobian[0], 0.5, 1e-10);  // dG/d[A] = 0.5
  EXPECT_NEAR(jacobian[1], 1.5, 1e-10);  // dG/d[B] = 1.5
}

TEST(LinearConstraint, JacobianIndependentOfConcentrations)
{
  // Test that Jacobian is constant (independent of concentrations) for linear constraints
  LinearConstraint constraint(
      "A_B_sum",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 2.0), StoichSpecies(Species("B"), 3.0) },
      1.0);

  // Test at two different concentration points
  std::vector<double> conc1 = { 0.1, 0.2 };
  std::vector<double> conc2 = { 10.0, 20.0 };
  std::vector<std::size_t> indices = { 0, 1 };

  std::vector<double> jacobian1(2);
  std::vector<double> jacobian2(2);

  constraint.Jacobian(conc1.data(), indices.data(), jacobian1.data());
  constraint.Jacobian(conc2.data(), indices.data(), jacobian2.data());

  // Jacobian should be the same regardless of concentrations
  EXPECT_NEAR(jacobian1[0], jacobian2[0], 1e-10);
  EXPECT_NEAR(jacobian1[1], jacobian2[1], 1e-10);
  EXPECT_NEAR(jacobian1[0], 2.0, 1e-10);
  EXPECT_NEAR(jacobian1[1], 3.0, 1e-10);
}
