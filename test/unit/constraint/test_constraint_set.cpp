// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/equilibrium_constraint.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <system_error>
#include <vector>

using namespace micm;

TEST(ConstraintSet, Construction)
{
  // Create a simple constraint set with one equilibrium constraint
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      1000.0));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 },
    { "B", 1 },
    { "AB", 2 }
  };

  std::size_t constraint_row_offset = 3;  // After the 3 species

  ConstraintSet set(std::move(constraints), variable_map, constraint_row_offset);

  EXPECT_EQ(set.Size(), 1);
}

TEST(ConstraintSet, NonZeroJacobianElements)
{
  // Create constraint set
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      1000.0));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 },
    { "B", 1 },
    { "AB", 2 }
  };

  std::size_t constraint_row_offset = 3;

  ConstraintSet set(std::move(constraints), variable_map, constraint_row_offset);

  auto non_zero_elements = set.NonZeroJacobianElements();

  // Constraint is at row 3 (after species A, B, AB)
  // It depends on A (col 0), B (col 1), AB (col 2)
  EXPECT_EQ(non_zero_elements.size(), 3);

  EXPECT_TRUE(non_zero_elements.count(std::make_pair(3, 0)));  // dG/dA at row 3, col 0
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(3, 1)));  // dG/dB at row 3, col 1
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(3, 2)));  // dG/dAB at row 3, col 2
}

TEST(ConstraintSet, MultipleConstraints)
{
  // Create constraint set with two equilibrium constraints
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      1000.0));
  constraints.push_back(EquilibriumConstraint(
      "C_D_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("C"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("D"), 1.0) },
      500.0));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 },
    { "B", 1 },
    { "AB", 2 },
    { "C", 3 },
    { "D", 4 }
  };

  std::size_t constraint_row_offset = 5;

  ConstraintSet set(std::move(constraints), variable_map, constraint_row_offset);

  EXPECT_EQ(set.Size(), 2);

  auto non_zero_elements = set.NonZeroJacobianElements();

  // First constraint at row 5: depends on A, B, AB
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(5, 0)));  // dG1/dA
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(5, 1)));  // dG1/dB
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(5, 2)));  // dG1/dAB

  // Second constraint at row 6: depends on C, D
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(6, 3)));  // dG2/dC
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(6, 4)));  // dG2/dD

  EXPECT_EQ(non_zero_elements.size(), 5);
}

TEST(ConstraintSet, AddForcingTerms)
{
  // Create constraint set
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      1000.0));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 },
    { "B", 1 },
    { "AB", 2 }
  };

  std::size_t num_species = 3;
  std::size_t num_constraints = 1;

  ConstraintSet set(std::move(constraints), variable_map, num_species);

  // State with 2 grid cells
  Matrix<double> state(2, num_species);
  state[0] = { 0.01, 0.01, 0.001 };  // Away from equilibrium
  state[1] = { 0.001, 0.001, 0.001 };  // At equilibrium

  // Extended forcing vector (species + constraints)
  Matrix<double> forcing(2, num_species + num_constraints, 0.0);

  set.AddForcingTerms(state, forcing);

  // For grid cell 0: G = 1000 * 0.01 * 0.01 - 0.001 = 0.1 - 0.001 = 0.099
  EXPECT_NEAR(forcing[0][3], 0.099, 1e-10);

  // For grid cell 1: G = 1000 * 0.001 * 0.001 - 0.001 = 0.001 - 0.001 = 0.0
  EXPECT_NEAR(forcing[1][3], 0.0, 1e-10);
}

TEST(ConstraintSet, SubtractJacobianTerms)
{
  // Create constraint set
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      1000.0));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 },
    { "B", 1 },
    { "AB", 2 }
  };

  std::size_t num_species = 3;
  std::size_t num_constraints = 1;
  std::size_t total_vars = num_species + num_constraints;

  ConstraintSet set(std::move(constraints), variable_map, num_species);

  // Get non-zero elements and build sparse Jacobian
  auto non_zero_elements = set.NonZeroJacobianElements();

  // Build a 4x4 sparse Jacobian (3 species + 1 constraint)
  // Include diagonal elements for species and constraint rows
  auto builder = SparseMatrix<double, SparseMatrixStandardOrdering>::Create(total_vars)
    .SetNumberOfBlocks(1)
    .InitialValue(0.0);

  for (std::size_t i = 0; i < total_vars; ++i)
    builder = builder.WithElement(i, i);  // Diagonals
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrix<double, SparseMatrixStandardOrdering> jacobian{ builder };

  set.SetJacobianFlatIds(jacobian);

  // State with 1 grid cell
  Matrix<double> state(1, num_species);
  state[0] = { 0.01, 0.02, 0.05 };

  set.SubtractJacobianTerms(state, jacobian);

  // G = K_eq * [A] * [B] - [AB]
  // dG/d[A] = K_eq * [B] = 1000 * 0.02 = 20
  // dG/d[B] = K_eq * [A] = 1000 * 0.01 = 10
  // dG/d[AB] = -1

  // Jacobian subtracts these values (matching ProcessSet convention)
  EXPECT_NEAR(jacobian[0][3][0], -20.0, 1e-10);  // J[constraint_row, A] -= dG/dA
  EXPECT_NEAR(jacobian[0][3][1], -10.0, 1e-10);  // J[constraint_row, B] -= dG/dB
  EXPECT_NEAR(jacobian[0][3][2], 1.0, 1e-10);    // J[constraint_row, AB] -= dG/dAB = -(-1) = 1
}

TEST(ConstraintSet, EmptyConstraintSet)
{
  // Empty constraint set should be valid and do nothing
  ConstraintSet set;

  EXPECT_EQ(set.Size(), 0);

  auto non_zero_elements = set.NonZeroJacobianElements();
  EXPECT_TRUE(non_zero_elements.empty());

  // AddForcingTerms and SubtractJacobianTerms should do nothing for empty set
  Matrix<double> state(1, 3);
  state[0] = { 0.1, 0.2, 0.3 };

  Matrix<double> forcing(1, 4, 1.0);

  set.AddForcingTerms(state, forcing);

  // Forcing should be unchanged
  EXPECT_DOUBLE_EQ(forcing[0][0], 1.0);
  EXPECT_DOUBLE_EQ(forcing[0][1], 1.0);
  EXPECT_DOUBLE_EQ(forcing[0][2], 1.0);
  EXPECT_DOUBLE_EQ(forcing[0][3], 1.0);
}

TEST(ConstraintSet, UnknownSpeciesThrows)
{
  // Creating a constraint with unknown species should throw
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "invalid",
      std::vector<StoichSpecies>{ StoichSpecies(Species("X"), 1.0), StoichSpecies(Species("Y"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("XY"), 1.0) },
      1000.0));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 },
    { "B", 1 }
  };

  EXPECT_THROW(
      ConstraintSet(std::move(constraints), variable_map, 2),
      std::system_error);
}

/// @brief Test 3D state (3 species) with 1 constraint
///
/// System: Species X, Y, Z with constraint X <-> Y (K_eq = 50)
/// State dimension: 3 species + 1 constraint = 4
/// Jacobian: 4x4 matrix
///
/// Constraint: G = K_eq * [X] - [Y] = 0
/// At equilibrium: [Y]/[X] = K_eq = 50
TEST(ConstraintSet, ThreeDStateOneConstraint)
{
  const double K_eq = 50.0;
  const std::size_t num_species = 3;
  const std::size_t num_constraints = 1;
  const std::size_t total_vars = num_species + num_constraints;

  // Create constraint: X <-> Y with K_eq = 50
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "X_Y_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("X"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("Y"), 1.0) },
      K_eq));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "X", 0 },
    { "Y", 1 },
    { "Z", 2 }
  };

  ConstraintSet set(std::move(constraints), variable_map, num_species);

  EXPECT_EQ(set.Size(), 1);

  // Check non-zero Jacobian elements
  auto non_zero_elements = set.NonZeroJacobianElements();
  EXPECT_EQ(non_zero_elements.size(), 2);  // dG/dX and dG/dY
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(3, 0)));  // dG/dX at row 3, col 0
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(3, 1)));  // dG/dY at row 3, col 1
  EXPECT_FALSE(non_zero_elements.count(std::make_pair(3, 2))); // Z not involved

  // Build sparse Jacobian (4x4)
  auto builder = SparseMatrix<double, SparseMatrixStandardOrdering>::Create(total_vars)
    .SetNumberOfBlocks(2)  // Test with 2 grid cells
    .InitialValue(0.0);

  for (std::size_t i = 0; i < total_vars; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrix<double, SparseMatrixStandardOrdering> jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  // State with 2 grid cells
  Matrix<double> state(2, num_species);
  // Grid cell 0: Away from equilibrium
  state[0][0] = 0.1;   // X
  state[0][1] = 2.0;   // Y
  state[0][2] = 0.5;   // Z (uninvolved)
  // Grid cell 1: At equilibrium (Y/X = 50)
  state[1][0] = 0.02;  // X
  state[1][1] = 1.0;   // Y = 50 * 0.02 = 1.0
  state[1][2] = 0.3;   // Z

  // Test forcing terms
  Matrix<double> forcing(2, total_vars, 0.0);
  set.AddForcingTerms(state, forcing);

  // Grid cell 0: G = K_eq * [X] - [Y] = 50 * 0.1 - 2.0 = 5.0 - 2.0 = 3.0
  EXPECT_NEAR(forcing[0][3], 3.0, 1e-10);
  // Grid cell 1: G = 50 * 0.02 - 1.0 = 1.0 - 1.0 = 0.0 (at equilibrium)
  EXPECT_NEAR(forcing[1][3], 0.0, 1e-10);

  // Test Jacobian terms
  set.SubtractJacobianTerms(state, jacobian);

  // For constraint G = K_eq * [X] - [Y]:
  // dG/dX = K_eq = 50
  // dG/dY = -1
  // Jacobian subtracts, so:
  // J[3,0] -= 50 => J[3,0] = -50
  // J[3,1] -= (-1) => J[3,1] = 1

  // Grid cell 0
  EXPECT_NEAR(jacobian[0][3][0], -K_eq, 1e-10);  // dG/dX
  EXPECT_NEAR(jacobian[0][3][1], 1.0, 1e-10);    // dG/dY (subtracted -1)
  // Grid cell 1
  EXPECT_NEAR(jacobian[1][3][0], -K_eq, 1e-10);
  EXPECT_NEAR(jacobian[1][3][1], 1.0, 1e-10);

  // Z column should be unaffected (diagonal only)
  EXPECT_NEAR(jacobian[0][2][2], 0.0, 1e-10);
  EXPECT_NEAR(jacobian[1][2][2], 0.0, 1e-10);
}

/// @brief Test 4D state (4 species) with 2 constraints
///
/// System: Species A, B, C, D with two constraints:
///   1. A <-> B with K_eq1 = 10
///   2. C + D <-> A with K_eq2 = 100
///
/// State dimension: 4 species + 2 constraints = 6
/// Jacobian: 6x6 matrix
///
/// Constraint 1: G1 = K_eq1 * [A] - [B] = 0
/// Constraint 2: G2 = K_eq2 * [C] * [D] - [A] = 0
TEST(ConstraintSet, FourDStateTwoConstraints)
{
  const double K_eq1 = 10.0;
  const double K_eq2 = 100.0;
  const std::size_t num_species = 4;
  const std::size_t num_constraints = 2;
  const std::size_t total_vars = num_species + num_constraints;

  // Create two constraints
  std::vector<Constraint> constraints;

  // Constraint 1: A <-> B with K_eq1 = 10
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
      K_eq1));

  // Constraint 2: C + D <-> A with K_eq2 = 100
  constraints.push_back(EquilibriumConstraint(
      "CD_A_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("C"), 1.0), StoichSpecies(Species("D"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
      K_eq2));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 },
    { "B", 1 },
    { "C", 2 },
    { "D", 3 }
  };

  ConstraintSet set(std::move(constraints), variable_map, num_species);

  EXPECT_EQ(set.Size(), 2);

  // Check non-zero Jacobian elements
  auto non_zero_elements = set.NonZeroJacobianElements();

  // Constraint 1 at row 4: depends on A (col 0), B (col 1)
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(4, 0)));  // dG1/dA
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(4, 1)));  // dG1/dB

  // Constraint 2 at row 5: depends on C (col 2), D (col 3), A (col 0)
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(5, 2)));  // dG2/dC
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(5, 3)));  // dG2/dD
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(5, 0)));  // dG2/dA

  // Total: 2 (from constraint 1) + 3 (from constraint 2) = 5
  EXPECT_EQ(non_zero_elements.size(), 5);

  // Build sparse Jacobian (6x6)
  auto builder = SparseMatrix<double, SparseMatrixStandardOrdering>::Create(total_vars)
    .SetNumberOfBlocks(3)  // Test with 3 grid cells
    .InitialValue(0.0);

  for (std::size_t i = 0; i < total_vars; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrix<double, SparseMatrixStandardOrdering> jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  // State with 3 grid cells
  Matrix<double> state(3, num_species);

  // Grid cell 0: Both constraints satisfied
  // If [A] = 0.1, then [B] = K_eq1 * [A] = 10 * 0.1 = 1.0
  // If [A] = 0.1 and K_eq2 * [C] * [D] = [A], then [C] * [D] = 0.1/100 = 0.001
  // Let [C] = 0.1, [D] = 0.01 => [C]*[D] = 0.001
  state[0][0] = 0.1;    // A
  state[0][1] = 1.0;    // B = 10 * 0.1
  state[0][2] = 0.1;    // C
  state[0][3] = 0.01;   // D, so C*D = 0.001, K_eq2*C*D = 0.1 = A

  // Grid cell 1: First constraint satisfied, second not
  state[1][0] = 0.2;    // A
  state[1][1] = 2.0;    // B = 10 * 0.2 (constraint 1 satisfied)
  state[1][2] = 0.1;    // C
  state[1][3] = 0.1;    // D, C*D = 0.01, K_eq2*C*D = 1.0 != 0.2

  // Grid cell 2: Neither constraint satisfied
  state[2][0] = 0.5;    // A
  state[2][1] = 3.0;    // B != 10 * 0.5 = 5.0
  state[2][2] = 0.2;    // C
  state[2][3] = 0.3;    // D, C*D = 0.06, K_eq2*C*D = 6.0 != 0.5

  // Test forcing terms
  Matrix<double> forcing(3, total_vars, 0.0);
  set.AddForcingTerms(state, forcing);

  // Grid cell 0: Both at equilibrium
  // G1 = K_eq1 * [A] - [B] = 10 * 0.1 - 1.0 = 0
  // G2 = K_eq2 * [C] * [D] - [A] = 100 * 0.1 * 0.01 - 0.1 = 0.1 - 0.1 = 0
  EXPECT_NEAR(forcing[0][4], 0.0, 1e-10);
  EXPECT_NEAR(forcing[0][5], 0.0, 1e-10);

  // Grid cell 1: First satisfied, second not
  // G1 = 10 * 0.2 - 2.0 = 0
  // G2 = 100 * 0.1 * 0.1 - 0.2 = 1.0 - 0.2 = 0.8
  EXPECT_NEAR(forcing[1][4], 0.0, 1e-10);
  EXPECT_NEAR(forcing[1][5], 0.8, 1e-10);

  // Grid cell 2: Neither satisfied
  // G1 = 10 * 0.5 - 3.0 = 5.0 - 3.0 = 2.0
  // G2 = 100 * 0.2 * 0.3 - 0.5 = 6.0 - 0.5 = 5.5
  EXPECT_NEAR(forcing[2][4], 2.0, 1e-10);
  EXPECT_NEAR(forcing[2][5], 5.5, 1e-10);

  // Test Jacobian terms
  set.SubtractJacobianTerms(state, jacobian);

  // Constraint 1: G1 = K_eq1 * [A] - [B]
  // dG1/dA = K_eq1 = 10
  // dG1/dB = -1

  // Constraint 2: G2 = K_eq2 * [C] * [D] - [A]
  // dG2/dC = K_eq2 * [D]
  // dG2/dD = K_eq2 * [C]
  // dG2/dA = -1

  // Grid cell 0:
  // Constraint 1 row (4):
  EXPECT_NEAR(jacobian[0][4][0], -K_eq1, 1e-10);  // -dG1/dA = -10
  EXPECT_NEAR(jacobian[0][4][1], 1.0, 1e-10);     // -dG1/dB = -(-1) = 1

  // Constraint 2 row (5):
  // dG2/dC = 100 * 0.01 = 1.0
  // dG2/dD = 100 * 0.1 = 10.0
  // dG2/dA = -1
  EXPECT_NEAR(jacobian[0][5][2], -K_eq2 * state[0][3], 1e-10);  // -dG2/dC
  EXPECT_NEAR(jacobian[0][5][3], -K_eq2 * state[0][2], 1e-10);  // -dG2/dD
  EXPECT_NEAR(jacobian[0][5][0], 1.0, 1e-10);                    // -dG2/dA = -(-1) = 1

  // Grid cell 1:
  EXPECT_NEAR(jacobian[1][4][0], -K_eq1, 1e-10);
  EXPECT_NEAR(jacobian[1][4][1], 1.0, 1e-10);
  EXPECT_NEAR(jacobian[1][5][2], -K_eq2 * state[1][3], 1e-10);
  EXPECT_NEAR(jacobian[1][5][3], -K_eq2 * state[1][2], 1e-10);
  EXPECT_NEAR(jacobian[1][5][0], 1.0, 1e-10);

  // Grid cell 2:
  EXPECT_NEAR(jacobian[2][4][0], -K_eq1, 1e-10);
  EXPECT_NEAR(jacobian[2][4][1], 1.0, 1e-10);
  EXPECT_NEAR(jacobian[2][5][2], -K_eq2 * state[2][3], 1e-10);
  EXPECT_NEAR(jacobian[2][5][3], -K_eq2 * state[2][2], 1e-10);
  EXPECT_NEAR(jacobian[2][5][0], 1.0, 1e-10);
}

/// @brief Test coupled constraints where constraints share species
///
/// System: A, B, C with constraints that both involve A:
///   1. A <-> B (K_eq1 = 5)
///   2. A <-> C (K_eq2 = 20)
///
/// This tests that the Jacobian correctly handles overlapping dependencies
TEST(ConstraintSet, CoupledConstraintsSharedSpecies)
{
  const double K_eq1 = 5.0;
  const double K_eq2 = 20.0;
  const std::size_t num_species = 3;
  const std::size_t num_constraints = 2;
  const std::size_t total_vars = num_species + num_constraints;

  std::vector<Constraint> constraints;

  // Both constraints depend on species A
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
      K_eq1));

  constraints.push_back(EquilibriumConstraint(
      "A_C_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("C"), 1.0) },
      K_eq2));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 },
    { "B", 1 },
    { "C", 2 }
  };

  ConstraintSet set(std::move(constraints), variable_map, num_species);

  EXPECT_EQ(set.Size(), 2);

  // Check Jacobian structure - both constraints depend on A
  auto non_zero_elements = set.NonZeroJacobianElements();

  // Constraint 1 (row 3): dG1/dA, dG1/dB
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(3, 0)));  // dG1/dA
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(3, 1)));  // dG1/dB

  // Constraint 2 (row 4): dG2/dA, dG2/dC
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(4, 0)));  // dG2/dA
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(4, 2)));  // dG2/dC

  EXPECT_EQ(non_zero_elements.size(), 4);

  // Build Jacobian
  auto builder = SparseMatrix<double, SparseMatrixStandardOrdering>::Create(total_vars)
    .SetNumberOfBlocks(1)
    .InitialValue(0.0);

  for (std::size_t i = 0; i < total_vars; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrix<double, SparseMatrixStandardOrdering> jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  // State at dual equilibrium: [B]/[A] = 5, [C]/[A] = 20
  // If [A] = 0.1, [B] = 0.5, [C] = 2.0
  Matrix<double> state(1, num_species);
  state[0][0] = 0.1;   // A
  state[0][1] = 0.5;   // B = 5 * 0.1
  state[0][2] = 2.0;   // C = 20 * 0.1

  // Test forcing terms
  Matrix<double> forcing(1, total_vars, 0.0);
  set.AddForcingTerms(state, forcing);

  // Both constraints should be satisfied
  // G1 = 5 * 0.1 - 0.5 = 0
  // G2 = 20 * 0.1 - 2.0 = 0
  EXPECT_NEAR(forcing[0][3], 0.0, 1e-10);
  EXPECT_NEAR(forcing[0][4], 0.0, 1e-10);

  // Test Jacobian terms
  set.SubtractJacobianTerms(state, jacobian);

  // Constraint 1: dG1/dA = 5, dG1/dB = -1
  EXPECT_NEAR(jacobian[0][3][0], -K_eq1, 1e-10);
  EXPECT_NEAR(jacobian[0][3][1], 1.0, 1e-10);

  // Constraint 2: dG2/dA = 20, dG2/dC = -1
  EXPECT_NEAR(jacobian[0][4][0], -K_eq2, 1e-10);
  EXPECT_NEAR(jacobian[0][4][2], 1.0, 1e-10);
}
