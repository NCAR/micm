// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/equilibrium_constraint.hpp>
#include <micm/system/species.hpp>
#include <micm/system/yield.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <system_error>
#include <vector>

using namespace micm;

TEST(EquilibriumConstraint, SimpleABEquilibrium)
{
  // Test: A + B <-> AB with K_eq = 1000
  // At equilibrium: [AB] / ([A][B]) = K_eq
  // Constraint: G = K_eq * [A] * [B] - [AB] = 0

  double K_eq = 1000.0;
  EquilibriumConstraint constraint(
      "A_B_equilibrium",
      std::vector<Yield>{ Yield(Species("A"), 1.0), Yield(Species("B"), 1.0) },  // reactants with stoich
      std::vector<Yield>{ Yield(Species("AB"), 1.0) },                           // products with stoich
      K_eq);

  EXPECT_EQ(constraint.name_, "A_B_equilibrium");
  EXPECT_EQ(constraint.species_dependencies_.size(), 3);
  EXPECT_EQ(constraint.species_dependencies_[0], "A");
  EXPECT_EQ(constraint.species_dependencies_[1], "B");
  EXPECT_EQ(constraint.species_dependencies_[2], "AB");

  // Test at equilibrium: [A] = 0.001, [B] = 0.001, [AB] = 0.001
  // K_eq * [A] * [B] = 1000 * 0.001 * 0.001 = 0.001 = [AB]
  // Residual should be 0
  std::vector<double> concentrations = { 0.001, 0.001, 0.001 };
  std::vector<std::size_t> indices = { 0, 1, 2 };

  double residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 0.0, 1e-10);

  // Test away from equilibrium: [A] = 0.01, [B] = 0.01, [AB] = 0.001
  // K_eq * [A] * [B] = 1000 * 0.01 * 0.01 = 0.1
  // Residual = 0.1 - 0.001 = 0.099
  concentrations = { 0.01, 0.01, 0.001 };
  residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 0.099, 1e-10);
}

TEST(EquilibriumConstraint, Jacobian)
{
  // Test Jacobian for A + B <-> AB with K_eq = 1000
  // G = K_eq * [A] * [B] - [AB]
  // dG/d[A] = K_eq * [B]
  // dG/d[B] = K_eq * [A]
  // dG/d[AB] = -1

  double K_eq = 1000.0;
  EquilibriumConstraint constraint(
      "A_B_equilibrium",
      std::vector<Yield>{ Yield(Species("A"), 1.0), Yield(Species("B"), 1.0) },
      std::vector<Yield>{ Yield(Species("AB"), 1.0) },
      K_eq);

  std::vector<double> concentrations = { 0.01, 0.02, 0.05 };
  std::vector<std::size_t> indices = { 0, 1, 2 };
  std::vector<double> jacobian(3);

  constraint.Jacobian(concentrations.data(), indices.data(), jacobian.data());

  EXPECT_NEAR(jacobian[0], K_eq * concentrations[1], 1e-10);  // dG/d[A] = K_eq * [B]
  EXPECT_NEAR(jacobian[1], K_eq * concentrations[0], 1e-10);  // dG/d[B] = K_eq * [A]
  EXPECT_NEAR(jacobian[2], -1.0, 1e-10);                      // dG/d[AB] = -1
}

TEST(EquilibriumConstraint, SingleReactantSingleProduct)
{
  // Test: A <-> B with K_eq = 10
  // At equilibrium: [B] / [A] = K_eq
  // Constraint: G = K_eq * [A] - [B] = 0

  double K_eq = 10.0;
  EquilibriumConstraint constraint(
      "A_B_simple",
      std::vector<Yield>{ Yield(Species("A"), 1.0) },
      std::vector<Yield>{ Yield(Species("B"), 1.0) },
      K_eq);

  EXPECT_EQ(constraint.species_dependencies_.size(), 2);

  // At equilibrium: [A] = 0.1, [B] = 1.0
  std::vector<double> concentrations = { 0.1, 1.0 };
  std::vector<std::size_t> indices = { 0, 1 };

  double residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 0.0, 1e-10);

  std::vector<double> jacobian(2);
  constraint.Jacobian(concentrations.data(), indices.data(), jacobian.data());
  EXPECT_NEAR(jacobian[0], K_eq, 1e-10);    // dG/d[A] = K_eq
  EXPECT_NEAR(jacobian[1], -1.0, 1e-10);    // dG/d[B] = -1
}

TEST(EquilibriumConstraint, TwoProductsOneReactant)
{
  // Test: 2A <-> B + C with K_eq = 100
  // At equilibrium: [B][C] / [A]^2 = K_eq
  // Constraint: G = K_eq * [A]^2 - [B] * [C] = 0

  double K_eq = 100.0;
  EquilibriumConstraint constraint(
      "dissociation",
      std::vector<Yield>{ Yield(Species("A"), 2.0) },                                // A with stoich 2
      std::vector<Yield>{ Yield(Species("B"), 1.0), Yield(Species("C"), 1.0) },
      K_eq);

  // Dependencies should be A, B, C
  EXPECT_EQ(constraint.species_dependencies_.size(), 3);

  // At equilibrium: [A] = 0.1, [B] = 0.5, [C] = 2.0
  // K_eq * [A]^2 = 100 * 0.01 = 1.0
  // [B] * [C] = 0.5 * 2.0 = 1.0
  std::vector<double> concentrations = { 0.1, 0.5, 2.0 };
  std::vector<std::size_t> indices = { 0, 1, 2 };

  double residual = constraint.Residual(concentrations.data(), indices.data());
  EXPECT_NEAR(residual, 0.0, 1e-10);
}

TEST(EquilibriumConstraint, InvalidEquilibriumConstant)
{
  // Test that negative or zero K_eq throws
  EXPECT_THROW(
      EquilibriumConstraint(
          "invalid",
          std::vector<Yield>{ Yield(Species("A"), 1.0) },
          std::vector<Yield>{ Yield(Species("B"), 1.0) },
          -1.0),
      std::system_error);

  EXPECT_THROW(
      EquilibriumConstraint(
          "invalid",
          std::vector<Yield>{ Yield(Species("A"), 1.0) },
          std::vector<Yield>{ Yield(Species("B"), 1.0) },
          0.0),
      std::system_error);
}
