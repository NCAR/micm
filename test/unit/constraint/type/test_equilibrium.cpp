// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>
#include <micm/util/jacobian_verification.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <system_error>
#include <vector>

using namespace micm;
using StandardSparseMatrix = SparseMatrix<double, SparseMatrixStandardOrdering>;

TEST(EquilibriumConstraint, Construction)
{
  // Test: A + B <-> AB with K_eq = 1000
  // At equilibrium: [AB] / ([A][B]) = K_eq
  // Constraint: G = K_eq * [A] * [B] - [AB] = 0

  double K_eq = 1000.0;
  EquilibriumConstraint constraint(
      "A_B_equilibrium",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      VantHoffParam{ .K_HLC_ref = K_eq, .delta_H = -2400.0 });

  EXPECT_EQ(constraint.name_, "A_B_equilibrium");
  EXPECT_EQ(constraint.species_dependencies_.size(), 3);
  EXPECT_EQ(constraint.species_dependencies_[0], "A");
  EXPECT_EQ(constraint.species_dependencies_[1], "B");
  EXPECT_EQ(constraint.species_dependencies_[2], "AB");
  EXPECT_EQ(constraint.reactants_.size(), 2);
  EXPECT_EQ(constraint.products_.size(), 1);
}

TEST(EquilibriumConstraint, AlgebraicSpecies)
{
  // Test that AlgebraicSpecies returns the first product species

  double K_eq = 1000.0;
  EquilibriumConstraint constraint(
      "A_B_equilibrium",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      VantHoffParam{ .K_HLC_ref = K_eq, .delta_H = -2400.0 });

  EXPECT_EQ(constraint.AlgebraicSpecies(), "AB");
}

TEST(EquilibriumConstraint, SingleReactantSingleProduct)
{
  // Test: A <-> B with K_eq = 10
  // At equilibrium: [B] / [A] = K_eq
  // Constraint: G = K_eq * [A] - [B] = 0

  double K_eq = 10.0;
  EquilibriumConstraint constraint(
      "A_B_simple",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
      VantHoffParam{ .K_HLC_ref = K_eq, .delta_H = -2400.0 });

  EXPECT_EQ(constraint.name_, "A_B_simple");
  EXPECT_EQ(constraint.species_dependencies_.size(), 2);
  EXPECT_EQ(constraint.species_dependencies_[0], "A");
  EXPECT_EQ(constraint.species_dependencies_[1], "B");
  EXPECT_EQ(constraint.AlgebraicSpecies(), "B");
}

TEST(EquilibriumConstraint, MultipleReactantsAndProducts)
{
  // Test: 2A <-> B + C with K_eq = 100
  // At equilibrium: [B][C] / [A]^2 = K_eq
  // Constraint: G = K_eq * [A]^2 - [B] * [C] = 0

  double K_eq = 100.0;
  EquilibriumConstraint constraint(
      "dissociation",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 2.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0), StoichSpecies(Species("C"), 1.0) },
      VantHoffParam{ .K_HLC_ref = K_eq, .delta_H = -2400.0 });

  EXPECT_EQ(constraint.name_, "dissociation");
  EXPECT_EQ(constraint.species_dependencies_.size(), 3);
  EXPECT_EQ(constraint.species_dependencies_[0], "A");
  EXPECT_EQ(constraint.species_dependencies_[1], "B");
  EXPECT_EQ(constraint.species_dependencies_[2], "C");
  EXPECT_EQ(constraint.reactants_.size(), 1);
  EXPECT_EQ(constraint.reactants_[0].coefficient_, 2.0);
  EXPECT_EQ(constraint.products_.size(), 2);
  EXPECT_EQ(constraint.AlgebraicSpecies(), "B");
}

TEST(EquilibriumConstraint, InvalidEquilibriumConstant)
{
  // Test that negative or zero K_eq throws
  EXPECT_THROW(
      EquilibriumConstraint(
          "invalid",
          std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
          std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
          VantHoffParam{ .K_HLC_ref = -1.0, .delta_H = -2400.0 }),
      micm::MicmException);

  EXPECT_THROW(
      EquilibriumConstraint(
          "invalid",
          std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
          std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
          VantHoffParam{ .K_HLC_ref = 0.0, .delta_H = -2400.0 }),
      micm::MicmException);
}

TEST(EquilibriumConstraint, EmptyReactantsThrows)
{
  EXPECT_THROW(
      EquilibriumConstraint(
          "invalid",
          std::vector<StoichSpecies>{},  // empty reactants
          std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
          VantHoffParam{ .K_HLC_ref = 1.0, .delta_H = -2400.0 }),
      micm::MicmException);
}

TEST(EquilibriumConstraint, EmptyProductsThrows)
{
  EXPECT_THROW(
      EquilibriumConstraint(
          "invalid",
          std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
          std::vector<StoichSpecies>{},  // empty products
          VantHoffParam{ .K_HLC_ref = 1.0, .delta_H = -2400.0 }),
      micm::MicmException);
}

TEST(EquilibriumConstraint, InvalidStoichiometryThrows)
{
  // Zero stoichiometry for reactant
  EXPECT_THROW(
      EquilibriumConstraint(
          "invalid",
          std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 0.0) },
          std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
          VantHoffParam{ .K_HLC_ref = 1.0, .delta_H = -2400.0 }),
      micm::MicmException);

  // Negative stoichiometry for reactant
  EXPECT_THROW(
      EquilibriumConstraint(
          "invalid",
          std::vector<StoichSpecies>{ StoichSpecies(Species("A"), -1.0) },
          std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
          VantHoffParam{ .K_HLC_ref = 1.0, .delta_H = -2400.0 }),
      micm::MicmException);

  // Zero stoichiometry for product
  EXPECT_THROW(
      EquilibriumConstraint(
          "invalid",
          std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
          std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 0.0) },
          VantHoffParam{ .K_HLC_ref = 1.0, .delta_H = -2400.0 }),
      micm::MicmException);

  // Negative stoichiometry for product
  EXPECT_THROW(
      EquilibriumConstraint(
          "invalid",
          std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
          std::vector<StoichSpecies>{ StoichSpecies(Species("B"), -2.0) },
          VantHoffParam{ .K_HLC_ref = 1.0, .delta_H = -2400.0 }),
      micm::MicmException);
}

// Integration tests using ConstraintSet to test residual and Jacobian computation

TEST(EquilibriumConstraint, ResidualComputationThroughConstraintSet)
{
  // Test: A + B <-> AB with K_eq = 1000
  // Constraint: G = K_eq * [A] * [B] - [AB] = 0

  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_equilibrium",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 1000.0, .delta_H = -2400.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "AB", 2 } };

  std::size_t num_species = 3;

  ConstraintSet<DenseMatrix, StandardSparseMatrix> set(std::move(constraints), variable_map);

  // Create sparse matrix for constraint setup
  auto non_zero_elements = set.NonZeroJacobianElements();

  auto builder = StandardSparseMatrix::Create(num_species).SetNumberOfBlocks(1).InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  StandardSparseMatrix jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);
  std::unordered_map<std::string, std::size_t> state_parameter_indices = { { "A_B_equilibrium", 0 } };
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  // Create state matrix with 1 grid cell and 3 species
  DenseMatrix state(1, 3);
  DenseMatrix forcing(1, 3);

  // Test at equilibrium: [A] = 0.001, [B] = 0.001, [AB] = 0.001
  // K_eq * [A] * [B] = 1000 * 0.001 * 0.001 = 0.001 = [AB]
  // Residual should be 0
  state[0][0] = 0.001;  // A
  state[0][1] = 0.001;  // B
  state[0][2] = 0.001;  // AB

  forcing.Fill(0.0);
  DenseMatrix state_parameters(1, 1, 1000.0);  // K_eq = 1000.0 for all grid cells
  set.AddForcingTerms(state, state_parameters, forcing);

  // The forcing term for AB (row 2) should be the constraint residual
  EXPECT_NEAR(forcing[0][2], 0.0, 1e-10);

  // Test away from equilibrium: [A] = 0.02, [B] = 0.03, [AB] = 0.05
  // K_eq * [A] * [B] = 1000 * 0.02 * 0.03 = 0.6
  // Residual = 0.6 - 0.05 = 0.55
  state[0][0] = 0.02;  // A
  state[0][1] = 0.03;  // B
  state[0][2] = 0.05;  // AB

  forcing.Fill(0.0);
  set.AddForcingTerms(state, state_parameters, forcing);

  EXPECT_NEAR(forcing[0][2], 0.55, 1e-10);
}

TEST(EquilibriumConstraint, JacobianComputationThroughConstraintSet)
{
  // Test Jacobian for A + B <-> AB with K_eq = 1000
  // G = K_eq * [A] * [B] - [AB]
  // dG/d[A] = K_eq * [B]
  // dG/d[B] = K_eq * [A]
  // dG/d[AB] = -1

  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_equilibrium",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 1000.0, .delta_H = -2400.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "AB", 2 } };

  std::size_t num_species = 3;

  ConstraintSet<DenseMatrix, StandardSparseMatrix> set(std::move(constraints), variable_map);

  // Create sparse matrix for Jacobian using builder
  auto non_zero_elements = set.NonZeroJacobianElements();

  auto builder = StandardSparseMatrix::Create(num_species).SetNumberOfBlocks(1).InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);  // Diagonals
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  StandardSparseMatrix jacobian{ builder };

  set.SetJacobianFlatIds(jacobian);
  std::unordered_map<std::string, std::size_t> state_parameter_indices = { { "A_B_equilibrium", 0 } };
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  // Create state matrix
  DenseMatrix state(1, 3);
  state[0][0] = 0.01;  // A
  state[0][1] = 0.02;  // B
  state[0][2] = 0.05;  // AB

  // Compute Jacobian
  DenseMatrix state_parameters(1, 1, 1000.0);  // K_eq = 1000.0
  set.SubtractJacobianTerms(state, state_parameters, jacobian);

  // The Jacobian computation uses subtraction convention
  // Row 2 (AB, the algebraic species): contains dG/d[A], dG/d[B], dG/d[AB]
  double dG_dA = 1000.0 * 0.02;  // K_eq * [B] = 20.0
  double dG_dB = 1000.0 * 0.01;  // K_eq * [A] = 10.0
  double dG_dAB = -1.0;

  // Due to SubtractJacobianTerms convention, the values are negated
  EXPECT_NEAR(jacobian[0][2][0], -dG_dA, 1e-8);   // -dG/d[A] = -20.0
  EXPECT_NEAR(jacobian[0][2][1], -dG_dB, 1e-8);   // -dG/d[B] = -10.0
  EXPECT_NEAR(jacobian[0][2][2], -dG_dAB, 1e-8);  // -dG/d[AB] = 1.0
}

TEST(EquilibriumConstraint, ComplexStoichiometryResidual)
{
  // Test: 2A <-> B + C with K_eq = 100
  // Constraint: G = K_eq * [A]^2 - [B] * [C] = 0

  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "dissociation",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 2.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0), StoichSpecies(Species("C"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 100.0, .delta_H = -2400.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "C", 2 } };

  std::size_t num_species = 3;

  ConstraintSet<DenseMatrix, StandardSparseMatrix> set(std::move(constraints), variable_map);

  // Create sparse matrix for constraint setup
  auto non_zero_elements = set.NonZeroJacobianElements();

  auto builder = StandardSparseMatrix::Create(num_species).SetNumberOfBlocks(1).InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  StandardSparseMatrix jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);
  std::unordered_map<std::string, std::size_t> state_parameter_indices = { { "dissociation", 0 } };
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  DenseMatrix state(1, 3);
  DenseMatrix forcing(1, 3);

  // At equilibrium: [A] = 0.1, [B] = 0.5, [C] = 2.0
  // K_eq * [A]^2 = 100 * 0.01 = 1.0
  // [B] * [C] = 0.5 * 2.0 = 1.0
  state[0][0] = 0.1;  // A
  state[0][1] = 0.5;  // B
  state[0][2] = 2.0;  // C

  forcing.Fill(0.0);
  DenseMatrix state_parameters(1, 1, 100.0);  // K_eq = 100.0
  set.AddForcingTerms(state, state_parameters, forcing);

  // The forcing term for B (row 1, first product) should be the constraint residual
  EXPECT_NEAR(forcing[0][1], 0.0, 1e-10);
}

TEST(EquilibriumConstraint, FiniteDifferenceJacobianSimple)
{
  // A + B <-> AB, K_eq = 1000 (at 298.15 K with delta_H = 0)
  // G = K_eq * [A] * [B] - [AB]
  // dG/dA = K_eq * [B], dG/dB = K_eq * [A], dG/dAB = -1
  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 1000.0, .delta_H = 0.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "AB", 2 } };
  const std::size_t num_species = 3;

  ConstraintSet<DenseMatrix, StandardSparseMatrix> set(std::move(constraints), variable_map);

  auto non_zero_elements = set.NonZeroJacobianElements();
  auto builder = StandardSparseMatrix::Create(num_species).SetNumberOfBlocks(2).InitialValue(0.0);
  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);
  StandardSparseMatrix jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);
  std::unordered_map<std::string, std::size_t> state_parameter_indices = { { "eq", 0 } };
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  DenseMatrix variables(2, num_species, 0.0);
  variables[0][0] = 0.02;   // A
  variables[0][1] = 0.03;   // B
  variables[0][2] = 0.05;   // AB
  variables[1][0] = 0.1;
  variables[1][1] = 0.2;
  variables[1][2] = 5.0;

  DenseMatrix state_parameters(2, 1, 1000.0);  // K_eq = 1000

  // Analytical Jacobian
  set.SubtractJacobianTerms(variables, state_parameters, jacobian);

  // FD Jacobian
  auto fd_wrapper = [&](const DenseMatrix& vars, DenseMatrix& forcing)
  { set.AddForcingTerms(vars, state_parameters, forcing); };

  auto fd_jac = FiniteDifferenceJacobian<DenseMatrix>(fd_wrapper, variables, num_species);

  auto comparison = CompareJacobianToFiniteDifference<DenseMatrix, StandardSparseMatrix>(
      jacobian, fd_jac, num_species, /*atol=*/1e-7, /*rtol=*/1e-7);

  EXPECT_TRUE(comparison.passed) << "Equilibrium constraint Jacobian mismatch: block=" << comparison.worst_block
                                 << " row=" << comparison.worst_row << " col=" << comparison.worst_col
                                 << " analytical=" << comparison.worst_analytical << " fd=" << comparison.worst_fd;

  auto sparsity = CheckJacobianSparsityCompleteness<DenseMatrix, StandardSparseMatrix>(
      jacobian, fd_jac, num_species);

  EXPECT_TRUE(sparsity.passed) << "Missing sparsity at block=" << sparsity.worst_block
                               << " row=" << sparsity.worst_row << " col=" << sparsity.worst_col
                               << " fd_value=" << sparsity.worst_fd;
}

TEST(EquilibriumConstraint, FiniteDifferenceJacobianComplexStoichiometry)
{
  // 2A <-> B + C, K_eq = 100
  // G = K_eq * [A]^2 - [B] * [C]
  // dG/dA = K_eq * 2 * [A], dG/dB = -[C], dG/dC = -[B]
  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "dissociation",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 2.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0), StoichSpecies(Species("C"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 100.0, .delta_H = 0.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "C", 2 } };
  const std::size_t num_species = 3;

  ConstraintSet<DenseMatrix, StandardSparseMatrix> set(std::move(constraints), variable_map);

  auto non_zero_elements = set.NonZeroJacobianElements();
  auto builder = StandardSparseMatrix::Create(num_species).SetNumberOfBlocks(1).InitialValue(0.0);
  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);
  StandardSparseMatrix jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);
  std::unordered_map<std::string, std::size_t> state_parameter_indices = { { "dissociation", 0 } };
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  DenseMatrix variables(1, num_species, 0.0);
  variables[0][0] = 0.15;  // A
  variables[0][1] = 0.8;   // B
  variables[0][2] = 1.5;   // C

  DenseMatrix state_parameters(1, 1, 100.0);

  set.SubtractJacobianTerms(variables, state_parameters, jacobian);

  auto fd_wrapper = [&](const DenseMatrix& vars, DenseMatrix& forcing)
  { set.AddForcingTerms(vars, state_parameters, forcing); };

  auto fd_jac = FiniteDifferenceJacobian<DenseMatrix>(fd_wrapper, variables, num_species);

  auto comparison = CompareJacobianToFiniteDifference<DenseMatrix, StandardSparseMatrix>(
      jacobian, fd_jac, num_species, /*atol=*/1e-7, /*rtol=*/1e-7);

  EXPECT_TRUE(comparison.passed) << "Complex stoichiometry Jacobian mismatch: block=" << comparison.worst_block
                                 << " row=" << comparison.worst_row << " col=" << comparison.worst_col
                                 << " analytical=" << comparison.worst_analytical << " fd=" << comparison.worst_fd;

  auto sparsity = CheckJacobianSparsityCompleteness<DenseMatrix, StandardSparseMatrix>(
      jacobian, fd_jac, num_species);

  EXPECT_TRUE(sparsity.passed) << "Missing sparsity at block=" << sparsity.worst_block
                               << " row=" << sparsity.worst_row << " col=" << sparsity.worst_col
                               << " fd_value=" << sparsity.worst_fd;
}
