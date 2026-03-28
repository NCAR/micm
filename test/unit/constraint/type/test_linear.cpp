// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/linear_constraint.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>
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

TEST(LinearConstraint, Construction)
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
}

TEST(LinearConstraint, AlgebraicSpecies)
{
  // Test that AlgebraicSpecies returns the last species in terms list
  LinearConstraint constraint(
      "A_B_C_conservation",
      std::vector<StoichSpecies>{
          StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0), StoichSpecies(Species("C"), 1.0) },
      10.0);

  EXPECT_EQ(constraint.AlgebraicSpecies(), "C");
}

TEST(LinearConstraint, WeightedTerms)
{
  // Test: 2*A + 3*B - C = 5.0
  // Constraint: G = 2*[A] + 3*[B] - [C] - 5.0 = 0

  LinearConstraint constraint(
      "weighted_sum",
      std::vector<StoichSpecies>{
          StoichSpecies(Species("A"), 2.0), StoichSpecies(Species("B"), 3.0), StoichSpecies(Species("C"), -1.0) },
      5.0);

  EXPECT_EQ(constraint.name_, "weighted_sum");
  EXPECT_EQ(constraint.species_dependencies_.size(), 3);
  EXPECT_EQ(constraint.constant_, 5.0);
  EXPECT_EQ(constraint.terms_.size(), 3);
  EXPECT_EQ(constraint.terms_[0].coefficient_, 2.0);
  EXPECT_EQ(constraint.terms_[1].coefficient_, 3.0);
  EXPECT_EQ(constraint.terms_[2].coefficient_, -1.0);
}

TEST(LinearConstraint, ZeroConstant)
{
  // Test: A - B = 0 (species balance)
  // Constraint: G = [A] - [B] = 0

  LinearConstraint constraint(
      "A_equals_B", std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), -1.0) }, 0.0);

  EXPECT_EQ(constraint.name_, "A_equals_B");
  EXPECT_EQ(constraint.constant_, 0.0);
  EXPECT_EQ(constraint.species_dependencies_.size(), 2);
}

TEST(LinearConstraint, FractionalCoefficients)
{
  // Test: 0.5*A + 1.5*B = 2.0
  // Constraint: G = 0.5*[A] + 1.5*[B] - 2.0 = 0

  LinearConstraint constraint(
      "fractional_conservation",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 0.5), StoichSpecies(Species("B"), 1.5) },
      2.0);

  EXPECT_EQ(constraint.name_, "fractional_conservation");
  EXPECT_EQ(constraint.constant_, 2.0);
  EXPECT_EQ(constraint.terms_[0].coefficient_, 0.5);
  EXPECT_EQ(constraint.terms_[1].coefficient_, 1.5);
}

// Integration tests using ConstraintSet to test residual and Jacobian computation

TEST(LinearConstraint, ResidualComputationThroughConstraintSet)
{
  // Test: A + B = 1.0 (total concentration conservation)
  // Constraint: G = [A] + [B] - 1.0 = 0

  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(LinearConstraint(
      "A_B_conservation",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      1.0));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 } };

  std::size_t num_species = 2;

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
  std::unordered_map<std::string, std::size_t> state_parameter_indices;  // Empty for linear constraints
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  // Create state matrix with 1 grid cell and 2 species
  DenseMatrix state(1, 2);
  DenseMatrix forcing(1, 2);

  // Test when constraint is satisfied: [A] = 0.3, [B] = 0.7
  // Residual = 0.3 + 0.7 - 1.0 = 0
  state[0][0] = 0.3;  // A
  state[0][1] = 0.7;  // B

  forcing.Fill(0.0);
  DenseMatrix state_parameters(1, 0);  // No parameters for linear constraints
  set.AddForcingTerms(state, state_parameters, forcing);

  // The forcing term for B (row 1, last term) should be the constraint residual
  EXPECT_NEAR(forcing[0][1], 0.0, 1e-10);

  // Test when constraint is not satisfied: [A] = 0.5, [B] = 0.6
  // Residual = 0.5 + 0.6 - 1.0 = 0.1
  state[0][0] = 0.5;  // A
  state[0][1] = 0.6;  // B

  forcing.Fill(0.0);
  set.AddForcingTerms(state, state_parameters, forcing);

  EXPECT_NEAR(forcing[0][1], 0.1, 1e-10);
}

TEST(LinearConstraint, JacobianComputationThroughConstraintSet)
{
  // Test Jacobian for A + B = 1.0
  // G = [A] + [B] - 1.0
  // dG/d[A] = 1.0
  // dG/d[B] = 1.0

  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(LinearConstraint(
      "A_B_conservation",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      1.0));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 } };

  std::size_t num_species = 2;

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
  std::unordered_map<std::string, std::size_t> state_parameter_indices;  // Empty for linear constraints
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  // Create state matrix
  DenseMatrix state(1, 2);
  state[0][0] = 0.3;  // A
  state[0][1] = 0.7;  // B

  // Compute Jacobian
  DenseMatrix state_parameters(1, 0);  // No parameters for linear constraints
  set.SubtractJacobianTerms(state, state_parameters, jacobian);

  // Due to SubtractJacobianTerms convention, the values are negated
  // Row 1 (B, the algebraic species): contains dG/d[A], dG/d[B]
  EXPECT_NEAR(jacobian[0][1][0], -1.0, 1e-10);  // -dG/d[A] = -1.0
  EXPECT_NEAR(jacobian[0][1][1], -1.0, 1e-10);  // -dG/d[B] = -1.0
}

TEST(LinearConstraint, WeightedSumResidualAndJacobian)
{
  // Test: 2*A + 3*B - C = 5.0
  // Constraint: G = 2*[A] + 3*[B] - [C] - 5.0 = 0

  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(LinearConstraint(
      "weighted_sum",
      std::vector<StoichSpecies>{
          StoichSpecies(Species("A"), 2.0), StoichSpecies(Species("B"), 3.0), StoichSpecies(Species("C"), -1.0) },
      5.0));

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
  std::unordered_map<std::string, std::size_t> state_parameter_indices;  // Empty for linear constraints
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  DenseMatrix state(1, 3);
  DenseMatrix forcing(1, 3);

  // Test when constraint is satisfied: [A] = 1.0, [B] = 2.0, [C] = 3.0
  // Residual = 2*1.0 + 3*2.0 - 3.0 - 5.0 = 2 + 6 - 3 - 5 = 0
  state[0][0] = 1.0;  // A
  state[0][1] = 2.0;  // B
  state[0][2] = 3.0;  // C

  forcing.Fill(0.0);
  DenseMatrix state_parameters(1, 0);  // No parameters for linear constraints
  set.AddForcingTerms(state, state_parameters, forcing);

  // The forcing term for C (row 2, last term) should be the constraint residual
  EXPECT_NEAR(forcing[0][2], 0.0, 1e-10);

  // Test Jacobian - reset jacobian values to zero first
  for (auto& elem : jacobian.AsVector())
  {
    elem = 0.0;
  }

  set.SubtractJacobianTerms(state, state_parameters, jacobian);

  // Row 2 (C, the algebraic species): contains dG/d[A], dG/d[B], dG/d[C]
  // Due to SubtractJacobianTerms convention, the values are negated
  EXPECT_NEAR(jacobian[0][2][0], -2.0, 1e-10);  // -dG/d[A] = -2.0
  EXPECT_NEAR(jacobian[0][2][1], -3.0, 1e-10);  // -dG/d[B] = -3.0
  EXPECT_NEAR(jacobian[0][2][2], 1.0, 1e-10);   // -dG/d[C] = -(-1.0) = 1.0
}

TEST(LinearConstraint, ThreeSpeciesConservationResidual)
{
  // Test: A + B + C = 10.0
  // Constraint: G = [A] + [B] + [C] - 10.0 = 0

  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(LinearConstraint(
      "ABC_total",
      std::vector<StoichSpecies>{
          StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0), StoichSpecies(Species("C"), 1.0) },
      10.0));

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
  std::unordered_map<std::string, std::size_t> state_parameter_indices;  // Empty for linear constraints
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  DenseMatrix state(1, 3);
  DenseMatrix forcing(1, 3);

  // Test when constraint is satisfied: [A] = 2.0, [B] = 3.0, [C] = 5.0
  state[0][0] = 2.0;  // A
  state[0][1] = 3.0;  // B
  state[0][2] = 5.0;  // C

  forcing.Fill(0.0);
  DenseMatrix state_parameters(1, 0);  // No parameters for linear constraints
  set.AddForcingTerms(state, state_parameters, forcing);

  EXPECT_NEAR(forcing[0][2], 0.0, 1e-10);

  // Test when constraint is violated: [A] = 2.0, [B] = 3.0, [C] = 6.0
  // Residual = 2.0 + 3.0 + 6.0 - 10.0 = 1.0
  state[0][0] = 2.0;  // A
  state[0][1] = 3.0;  // B
  state[0][2] = 6.0;  // C

  forcing.Fill(0.0);
  set.AddForcingTerms(state, state_parameters, forcing);

  EXPECT_NEAR(forcing[0][2], 1.0, 1e-10);
}

TEST(LinearConstraint, ZeroConstantResidual)
{
  // Test: A - B = 0 (species balance)
  // Constraint: G = [A] - [B] = 0

  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(LinearConstraint(
      "A_equals_B", std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), -1.0) }, 0.0));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 } };

  std::size_t num_species = 2;

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
  std::unordered_map<std::string, std::size_t> state_parameter_indices;  // Empty for linear constraints
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  DenseMatrix state(1, 2);
  DenseMatrix forcing(1, 2);

  // Test when constraint is satisfied: [A] = 0.5, [B] = 0.5
  state[0][0] = 0.5;  // A
  state[0][1] = 0.5;  // B

  forcing.Fill(0.0);
  DenseMatrix state_parameters(1, 0);  // No parameters for linear constraints
  set.AddForcingTerms(state, state_parameters, forcing);

  EXPECT_NEAR(forcing[0][1], 0.0, 1e-10);

  // Test when constraint is violated: [A] = 0.6, [B] = 0.5
  // Residual = 0.6 - 0.5 = 0.1
  state[0][0] = 0.6;  // A
  state[0][1] = 0.5;  // B

  forcing.Fill(0.0);
  set.AddForcingTerms(state, state_parameters, forcing);

  EXPECT_NEAR(forcing[0][1], 0.1, 1e-10);
}

TEST(LinearConstraint, FractionalCoefficientsResidualAndJacobian)
{
  // Test: 0.5*A + 1.5*B = 2.0
  // Constraint: G = 0.5*[A] + 1.5*[B] - 2.0 = 0

  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(LinearConstraint(
      "fractional_conservation",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 0.5), StoichSpecies(Species("B"), 1.5) },
      2.0));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 } };

  std::size_t num_species = 2;

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
  std::unordered_map<std::string, std::size_t> state_parameter_indices;  // Empty for linear constraints
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  DenseMatrix state(1, 2);
  DenseMatrix forcing(1, 2);

  // Test when constraint is satisfied: [A] = 1.0, [B] = 1.0
  // 0.5*1.0 + 1.5*1.0 = 0.5 + 1.5 = 2.0 ✓
  state[0][0] = 1.0;  // A
  state[0][1] = 1.0;  // B

  forcing.Fill(0.0);
  DenseMatrix state_parameters(1, 0);  // No parameters for linear constraints
  set.AddForcingTerms(state, state_parameters, forcing);

  EXPECT_NEAR(forcing[0][1], 0.0, 1e-10);

  // Test Jacobian computation
  set.SubtractJacobianTerms(state, state_parameters, jacobian);

  // Due to SubtractJacobianTerms convention, the values are negated
  EXPECT_NEAR(jacobian[0][1][0], -0.5, 1e-10);  // -dG/d[A] = -0.5
  EXPECT_NEAR(jacobian[0][1][1], -1.5, 1e-10);  // -dG/d[B] = -1.5
}

TEST(LinearConstraint, JacobianIndependentOfConcentrations)
{
  // Test that Jacobian is constant (independent of concentrations) for linear constraints

  using DenseMatrix = Matrix<double>;

  std::vector<Constraint> constraints;
  constraints.push_back(LinearConstraint(
      "A_B_sum", std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 2.0), StoichSpecies(Species("B"), 3.0) }, 1.0));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 } };

  std::size_t num_species = 2;

  ConstraintSet<DenseMatrix, StandardSparseMatrix> set(std::move(constraints), variable_map);

  auto non_zero_elements = set.NonZeroJacobianElements();

  auto builder = StandardSparseMatrix::Create(num_species).SetNumberOfBlocks(1).InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  StandardSparseMatrix jacobian{ builder };

  set.SetJacobianFlatIds(jacobian);
  std::unordered_map<std::string, std::size_t> state_parameter_indices;  // Empty for linear constraints
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  // Test at two different concentration points
  DenseMatrix state1(1, 2);
  state1[0][0] = 0.1;  // A
  state1[0][1] = 0.2;  // B

  DenseMatrix state2(1, 2);
  state2[0][0] = 10.0;  // A
  state2[0][1] = 20.0;  // B

  // Compute Jacobian at first state
  DenseMatrix state_parameters(1, 0);  // No parameters for linear constraints
  set.SubtractJacobianTerms(state1, state_parameters, jacobian);

  double jac1_dA = jacobian[0][1][0];
  double jac1_dB = jacobian[0][1][1];

  // Reset jacobian and compute at second state
  for (auto& elem : jacobian.AsVector())
  {
    elem = 0.0;
  }

  set.SubtractJacobianTerms(state2, state_parameters, jacobian);

  double jac2_dA = jacobian[0][1][0];
  double jac2_dB = jacobian[0][1][1];

  // Jacobian should be the same regardless of concentrations
  EXPECT_NEAR(jac1_dA, jac2_dA, 1e-10);
  EXPECT_NEAR(jac1_dB, jac2_dB, 1e-10);
  EXPECT_NEAR(jac1_dA, -2.0, 1e-10);  // -dG/d[A] = -2.0
  EXPECT_NEAR(jac1_dB, -3.0, 1e-10);  // -dG/d[B] = -3.0
}