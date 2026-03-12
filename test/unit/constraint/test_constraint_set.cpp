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
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

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

  ConstraintSet set(std::move(constraints), variable_map);

  EXPECT_EQ(set.Size(), 1);
}

TEST(ConstraintSet, ReplaceStateRowsMapsToAlgebraicSpecies)
{
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "B_C_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("C"), 1.0) },
      10.0));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 },
    { "B", 1 },
    { "C", 2 }
  };

  ConstraintSet set(std::move(constraints), variable_map);

  EXPECT_EQ(set.Size(), 1);
  EXPECT_EQ(set.AlgebraicVariableIds().size(), 1);
  EXPECT_TRUE(set.AlgebraicVariableIds().count(2));  // C row is algebraic

  auto non_zero_elements = set.NonZeroJacobianElements();
  EXPECT_EQ(non_zero_elements.size(), 2);
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 1)));  // row C, col B
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 2)));  // row C, col C
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

  ConstraintSet set(std::move(constraints), variable_map);

  auto non_zero_elements = set.NonZeroJacobianElements();

  // Algebraic species = AB (index 2), constraint replaces row 2
  // Dependencies: A (col 0), B (col 1), AB (col 2) + diagonal (2,2)
  EXPECT_EQ(non_zero_elements.size(), 3);

  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 0)));  // dG/dA at row 2, col 0
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 1)));  // dG/dB at row 2, col 1
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 2)));  // dG/dAB at row 2, col 2 (+ diagonal)
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

  ConstraintSet set(std::move(constraints), variable_map);

  EXPECT_EQ(set.Size(), 2);

  auto non_zero_elements = set.NonZeroJacobianElements();

  // First constraint: algebraic = AB (index 2), replaces row 2
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 0)));  // dG1/dA
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 1)));  // dG1/dB
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 2)));  // dG1/dAB + diagonal

  // Second constraint: algebraic = D (index 4), replaces row 4
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(4, 3)));  // dG2/dC
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(4, 4)));  // dG2/dD + diagonal

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

  ConstraintSet set(std::move(constraints), variable_map);

  // State with 2 grid cells
  Matrix<double> state(2, num_species);
  state[0] = { 0.01, 0.01, 0.001 };  // Away from equilibrium
  state[1] = { 0.001, 0.001, 0.001 };  // At equilibrium

  // Forcing vector (same size as state; constraint replaces AB row at index 2)
  Matrix<double> forcing(2, num_species, 0.0);

  set.AddForcingTerms(state, forcing);

  // For grid cell 0: G = 1000 * 0.01 * 0.01 - 0.001 = 0.1 - 0.001 = 0.099
  EXPECT_NEAR(forcing[0][2], 0.099, 1e-10);

  // For grid cell 1: G = 1000 * 0.001 * 0.001 - 0.001 = 0.001 - 0.001 = 0.0
  EXPECT_NEAR(forcing[1][2], 0.0, 1e-10);
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

  ConstraintSet set(std::move(constraints), variable_map);

  // Get non-zero elements and build sparse Jacobian
  auto non_zero_elements = set.NonZeroJacobianElements();

  // Build a 3x3 sparse Jacobian (constraint replaces AB's row)
  auto builder = SparseMatrix<double, SparseMatrixStandardOrdering>::Create(num_species)
    .SetNumberOfBlocks(1)
    .InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
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
  // Constraint replaces row 2 (AB's row)
  EXPECT_NEAR(jacobian[0][2][0], -20.0, 1e-10);  // J[2, A] -= dG/dA
  EXPECT_NEAR(jacobian[0][2][1], -10.0, 1e-10);  // J[2, B] -= dG/dB
  EXPECT_NEAR(jacobian[0][2][2], 1.0, 1e-10);    // J[2, AB] -= dG/dAB = -(-1) = 1
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
      ConstraintSet(std::move(constraints), variable_map),
      std::system_error);
}

/// @brief Test 3D state (3 species) with 1 constraint
///
/// System: Species X, Y, Z with constraint X <-> Y (K_eq = 50)
/// Algebraic species = Y (first product), replaces row 1
/// Jacobian: 3x3 matrix
///
/// Constraint: G = K_eq * [X] - [Y] = 0
/// At equilibrium: [Y]/[X] = K_eq = 50
TEST(ConstraintSet, ThreeDStateOneConstraint)
{
  const double K_eq = 50.0;
  const std::size_t num_species = 3;

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

  ConstraintSet set(std::move(constraints), variable_map);

  EXPECT_EQ(set.Size(), 1);

  // Check non-zero Jacobian elements
  // Algebraic species = Y (index 1), constraint replaces row 1
  auto non_zero_elements = set.NonZeroJacobianElements();
  EXPECT_EQ(non_zero_elements.size(), 2);  // (1,0) dG/dX, (1,1) dG/dY + diagonal
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(1, 0)));  // dG/dX at row 1, col 0
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(1, 1)));  // dG/dY at row 1, col 1 + diagonal
  EXPECT_FALSE(non_zero_elements.count(std::make_pair(1, 2))); // Z not involved

  // Build sparse Jacobian (3x3)
  auto builder = SparseMatrix<double, SparseMatrixStandardOrdering>::Create(num_species)
    .SetNumberOfBlocks(2)  // Test with 2 grid cells
    .InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
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
  Matrix<double> forcing(2, num_species, 0.0);
  set.AddForcingTerms(state, forcing);

  // Constraint replaces row 1 (Y's row)
  // Grid cell 0: G = K_eq * [X] - [Y] = 50 * 0.1 - 2.0 = 3.0
  EXPECT_NEAR(forcing[0][1], 3.0, 1e-10);
  // Grid cell 1: G = 50 * 0.02 - 1.0 = 0.0 (at equilibrium)
  EXPECT_NEAR(forcing[1][1], 0.0, 1e-10);

  // Test Jacobian terms
  set.SubtractJacobianTerms(state, jacobian);

  // For constraint G = K_eq * [X] - [Y]:
  // dG/dX = K_eq = 50
  // dG/dY = -1
  // Jacobian subtracts at row 1:

  // Grid cell 0
  EXPECT_NEAR(jacobian[0][1][0], -K_eq, 1e-10);  // dG/dX
  EXPECT_NEAR(jacobian[0][1][1], 1.0, 1e-10);    // dG/dY (subtracted -1)
  // Grid cell 1
  EXPECT_NEAR(jacobian[1][1][0], -K_eq, 1e-10);
  EXPECT_NEAR(jacobian[1][1][1], 1.0, 1e-10);

  // Z row should be unaffected
  EXPECT_NEAR(jacobian[0][2][2], 0.0, 1e-10);
  EXPECT_NEAR(jacobian[1][2][2], 0.0, 1e-10);
}

/// @brief Test 4D state (4 species) with 2 constraints
///
/// System: Species A, B, C, D with two constraints:
///   1. A <-> B with K_eq1 = 10, algebraic species = B (row 1)
///   2. C + D <-> A with K_eq2 = 100, algebraic species = A (row 0)
/// Jacobian: 4x4 matrix
///
/// Constraint 1: G1 = K_eq1 * [A] - [B] = 0
/// Constraint 2: G2 = K_eq2 * [C] * [D] - [A] = 0
TEST(ConstraintSet, FourDStateTwoConstraints)
{
  const double K_eq1 = 10.0;
  const double K_eq2 = 100.0;
  const std::size_t num_species = 4;

  // Create two constraints
  std::vector<Constraint> constraints;

  // Constraint 1: A <-> B with K_eq1 = 10, algebraic species = B (row 1)
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
      K_eq1));

  // Constraint 2: C + D <-> A with K_eq2 = 100, algebraic species = A (row 0)
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

  ConstraintSet set(std::move(constraints), variable_map);

  EXPECT_EQ(set.Size(), 2);

  // Check non-zero Jacobian elements
  auto non_zero_elements = set.NonZeroJacobianElements();

  // Constraint 1 replaces row 1 (B): depends on A (col 0), B (col 1)
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(1, 0)));  // dG1/dA
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(1, 1)));  // dG1/dB + diagonal

  // Constraint 2 replaces row 0 (A): depends on C (col 2), D (col 3), A (col 0)
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(0, 2)));  // dG2/dC
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(0, 3)));  // dG2/dD
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(0, 0)));  // dG2/dA + diagonal

  // Total: 2 (from constraint 1) + 3 (from constraint 2) = 5
  EXPECT_EQ(non_zero_elements.size(), 5);

  // Build sparse Jacobian (4x4)
  auto builder = SparseMatrix<double, SparseMatrixStandardOrdering>::Create(num_species)
    .SetNumberOfBlocks(3)  // Test with 3 grid cells
    .InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrix<double, SparseMatrixStandardOrdering> jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  // State with 3 grid cells
  Matrix<double> state(3, num_species);

  // Grid cell 0: Both constraints satisfied
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
  Matrix<double> forcing(3, num_species, 0.0);
  set.AddForcingTerms(state, forcing);

  // Constraint 1 replaces row 1, Constraint 2 replaces row 0
  // Grid cell 0: Both at equilibrium
  EXPECT_NEAR(forcing[0][1], 0.0, 1e-10);   // G1 = 10 * 0.1 - 1.0 = 0
  EXPECT_NEAR(forcing[0][0], 0.0, 1e-10);   // G2 = 100 * 0.1 * 0.01 - 0.1 = 0

  // Grid cell 1: First satisfied, second not
  EXPECT_NEAR(forcing[1][1], 0.0, 1e-10);   // G1 = 10 * 0.2 - 2.0 = 0
  EXPECT_NEAR(forcing[1][0], 0.8, 1e-10);   // G2 = 100 * 0.1 * 0.1 - 0.2 = 0.8

  // Grid cell 2: Neither satisfied
  EXPECT_NEAR(forcing[2][1], 2.0, 1e-10);   // G1 = 10 * 0.5 - 3.0 = 2.0
  EXPECT_NEAR(forcing[2][0], 5.5, 1e-10);   // G2 = 100 * 0.2 * 0.3 - 0.5 = 5.5

  // Test Jacobian terms
  set.SubtractJacobianTerms(state, jacobian);

  // Constraint 1 at row 1: dG1/dA = K_eq1, dG1/dB = -1
  // Constraint 2 at row 0: dG2/dC = K_eq2*[D], dG2/dD = K_eq2*[C], dG2/dA = -1

  // Grid cell 0:
  EXPECT_NEAR(jacobian[0][1][0], -K_eq1, 1e-10);
  EXPECT_NEAR(jacobian[0][1][1], 1.0, 1e-10);
  EXPECT_NEAR(jacobian[0][0][2], -K_eq2 * state[0][3], 1e-10);
  EXPECT_NEAR(jacobian[0][0][3], -K_eq2 * state[0][2], 1e-10);
  EXPECT_NEAR(jacobian[0][0][0], 1.0, 1e-10);

  // Grid cell 1:
  EXPECT_NEAR(jacobian[1][1][0], -K_eq1, 1e-10);
  EXPECT_NEAR(jacobian[1][1][1], 1.0, 1e-10);
  EXPECT_NEAR(jacobian[1][0][2], -K_eq2 * state[1][3], 1e-10);
  EXPECT_NEAR(jacobian[1][0][3], -K_eq2 * state[1][2], 1e-10);
  EXPECT_NEAR(jacobian[1][0][0], 1.0, 1e-10);

  // Grid cell 2:
  EXPECT_NEAR(jacobian[2][1][0], -K_eq1, 1e-10);
  EXPECT_NEAR(jacobian[2][1][1], 1.0, 1e-10);
  EXPECT_NEAR(jacobian[2][0][2], -K_eq2 * state[2][3], 1e-10);
  EXPECT_NEAR(jacobian[2][0][3], -K_eq2 * state[2][2], 1e-10);
  EXPECT_NEAR(jacobian[2][0][0], 1.0, 1e-10);
}

/// @brief Test coupled constraints where constraints share species
///
/// System: A, B, C with constraints that both involve A:
///   1. A <-> B (K_eq1 = 5), algebraic = B (row 1)
///   2. A <-> C (K_eq2 = 20), algebraic = C (row 2)
///
/// This tests that the Jacobian correctly handles overlapping dependencies
TEST(ConstraintSet, CoupledConstraintsSharedSpecies)
{
  const double K_eq1 = 5.0;
  const double K_eq2 = 20.0;
  const std::size_t num_species = 3;

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

  ConstraintSet set(std::move(constraints), variable_map);

  EXPECT_EQ(set.Size(), 2);

  // Check Jacobian structure - both constraints depend on A
  auto non_zero_elements = set.NonZeroJacobianElements();

  // Constraint 1 replaces row 1 (B): dG1/dA, dG1/dB + diagonal
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(1, 0)));  // dG1/dA
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(1, 1)));  // dG1/dB + diagonal

  // Constraint 2 replaces row 2 (C): dG2/dA, dG2/dC + diagonal
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 0)));  // dG2/dA
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 2)));  // dG2/dC + diagonal

  EXPECT_EQ(non_zero_elements.size(), 4);

  // Build Jacobian (3x3)
  auto builder = SparseMatrix<double, SparseMatrixStandardOrdering>::Create(num_species)
    .SetNumberOfBlocks(1)
    .InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrix<double, SparseMatrixStandardOrdering> jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  // State at dual equilibrium: [B]/[A] = 5, [C]/[A] = 20
  Matrix<double> state(1, num_species);
  state[0][0] = 0.1;   // A
  state[0][1] = 0.5;   // B = 5 * 0.1
  state[0][2] = 2.0;   // C = 20 * 0.1

  // Test forcing terms
  Matrix<double> forcing(1, num_species, 0.0);
  set.AddForcingTerms(state, forcing);

  // Both constraints should be satisfied
  EXPECT_NEAR(forcing[0][1], 0.0, 1e-10);  // G1 at row 1
  EXPECT_NEAR(forcing[0][2], 0.0, 1e-10);  // G2 at row 2

  // Test Jacobian terms
  set.SubtractJacobianTerms(state, jacobian);

  // Constraint 1 at row 1: dG1/dA = 5, dG1/dB = -1
  EXPECT_NEAR(jacobian[0][1][0], -K_eq1, 1e-10);
  EXPECT_NEAR(jacobian[0][1][1], 1.0, 1e-10);

  // Constraint 2 at row 2: dG2/dA = 20, dG2/dC = -1
  EXPECT_NEAR(jacobian[0][2][0], -K_eq2, 1e-10);
  EXPECT_NEAR(jacobian[0][2][2], 1.0, 1e-10);
}

TEST(ConstraintSet, VectorizedMatricesRespectGridCellIndexing)
{
  const std::size_t num_species = 3;

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

  ConstraintSet set(std::move(constraints), variable_map);
  auto non_zero_elements = set.NonZeroJacobianElements();

  // Constraint replaces AB's row (index 2), Jacobian is 3x3
  auto builder = SparseMatrix<double, SparseMatrixVectorOrdering<4>>::Create(num_species)
    .SetNumberOfBlocks(3)
    .InitialValue(0.0);
  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (const auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrix<double, SparseMatrixVectorOrdering<4>> jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  VectorMatrix<double, 4> state(3, num_species, 0.0);
  state[0] = { 0.01, 0.02, 0.05 };
  state[1] = { 0.03, 0.01, 0.2 };
  state[2] = { 0.001, 0.002, 0.004 };

  VectorMatrix<double, 4> forcing(3, num_species, 0.0);
  set.AddForcingTerms(state, forcing);

  // Constraint residual replaces row 2 (AB)
  EXPECT_NEAR(forcing[0][2], 1000.0 * 0.01 * 0.02 - 0.05, 1e-12);
  EXPECT_NEAR(forcing[1][2], 1000.0 * 0.03 * 0.01 - 0.2, 1e-12);
  EXPECT_NEAR(forcing[2][2], 1000.0 * 0.001 * 0.002 - 0.004, 1e-12);

  set.SubtractJacobianTerms(state, jacobian);

  // Jacobian entries at row 2 (AB's row, replaced by constraint)
  EXPECT_NEAR(jacobian[0][2][0], -(1000.0 * 0.02), 1e-12);
  EXPECT_NEAR(jacobian[0][2][1], -(1000.0 * 0.01), 1e-12);
  EXPECT_NEAR(jacobian[0][2][2], 1.0, 1e-12);

  EXPECT_NEAR(jacobian[1][2][0], -(1000.0 * 0.01), 1e-12);
  EXPECT_NEAR(jacobian[1][2][1], -(1000.0 * 0.03), 1e-12);
  EXPECT_NEAR(jacobian[1][2][2], 1.0, 1e-12);

  EXPECT_NEAR(jacobian[2][2][0], -(1000.0 * 0.002), 1e-12);
  EXPECT_NEAR(jacobian[2][2][1], -(1000.0 * 0.001), 1e-12);
  EXPECT_NEAR(jacobian[2][2][2], 1.0, 1e-12);
}
