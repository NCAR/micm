// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <micm/constraint/constraint.hpp>
#include <micm/constraint/constraint_set.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>
#include <micm/system/species.hpp>
#include <micm/system/stoich_species.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <set>
#include <system_error>
#include <unordered_map>
#include <vector>

using namespace micm;

template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testConstruction()
{
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "AB", 2 } };

  ConstraintSetPolicy set(std::move(constraints), variable_map);

  EXPECT_EQ(set.Size(), 1);
}

template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testReplaceStateRowsMapsToAlgebraicSpecies()
{
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "B_C_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("C"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "C", 2 } };

  ConstraintSetPolicy set(std::move(constraints), variable_map);

  EXPECT_EQ(set.Size(), 1);
  EXPECT_EQ(set.AlgebraicVariableIds().size(), 1);
  EXPECT_TRUE(set.AlgebraicVariableIds().count(2));  // C row is algebraic

  auto non_zero_elements = set.NonZeroJacobianElements();
  EXPECT_EQ(non_zero_elements.size(), 2);
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 1)));  // row C, col B
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 2)));  // row C, col C
}

template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testNonZeroJacobianElements()
{
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "AB", 2 } };

  ConstraintSetPolicy set(std::move(constraints), variable_map);

  auto non_zero_elements = set.NonZeroJacobianElements();

  // Algebraic species = AB (index 2), constraint replaces row 2
  // Dependencies: A (col 0), B (col 1), AB (col 2) + diagonal (2,2)
  EXPECT_EQ(non_zero_elements.size(), 3);
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 0)));  // dG/dA at row 2, col 0
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 1)));  // dG/dB at row 2, col 1
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(2, 2)));  // dG/dAB at row 2, col 2 (+ diagonal)
}

template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testMultipleConstraints()
{
  // Create constraint set with two equilibrium constraints
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));
  constraints.push_back(EquilibriumConstraint(
      "C_D_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("C"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("D"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = {
    { "A", 0 }, { "B", 1 }, { "AB", 2 }, { "C", 3 }, { "D", 4 }
  };

  ConstraintSetPolicy set(std::move(constraints), variable_map);

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

template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testAddForcingTerms()
{
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "AB", 2 } };

  std::size_t num_species = 3;

  ConstraintSetPolicy set(std::move(constraints), variable_map);

  // Build sparse Jacobian for SetConstraintFunctions
  auto non_zero_elements = set.NonZeroJacobianElements();
  auto builder = SparseMatrixPolicy::Create(num_species).SetNumberOfBlocks(2).InitialValue(0.0);
  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);
  SparseMatrixPolicy jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  std::unordered_map<std::string, std::size_t> state_parameter_indices = { { "A_B_eq", 0 } };
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  // State with 2 grid cells
  DenseMatrixPolicy state(2, num_species);
  state[0] = { 0.2, 0.4, 0.6 };  // Away from equilibrium
  state[1] = { 0.1, 0.3, 0.7 };  // Away from equilibrium

  // Forcing vector (same size as state; constraint replaces AB row at index 2)
  DenseMatrixPolicy forcing(2, num_species, 0.0);

  // State parameters: K_eq for each grid cell (2 cells, 1 parameter)
  DenseMatrixPolicy state_parameters(2, 1, 3.3e-2);

  set.AddForcingTerms(state, state_parameters, forcing);

  // For grid cell 0: G = K_eq * 0.2 * 0.4 - 0.6 = 3.3e-2 * 0.08 - 0.6 = 0.00264 - 0.6 = -0.59736
  EXPECT_NEAR(forcing[0][2], -0.59736, 1e-5);

  // For grid cell 1: G = K_eq * 0.1 * 0.3 - 0.7 = 3.3e-2 * 0.03 - 0.7 = 0.00099 - 0.7 = -0.69901
  EXPECT_NEAR(forcing[1][2], -0.69901, 1e-5);
}

template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testSubtractJacobianTerms()
{
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "AB", 2 } };

  std::size_t num_species = 3;

  ConstraintSetPolicy set(std::move(constraints), variable_map);

  // Get non-zero elements and build sparse Jacobian
  auto non_zero_elements = set.NonZeroJacobianElements();

  // Build a 3x3 sparse Jacobian (constraint replaces AB's row)
  auto builder = SparseMatrixPolicy::Create(num_species).SetNumberOfBlocks(1).InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);  // Diagonals
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrixPolicy jacobian{ builder };

  set.SetJacobianFlatIds(jacobian);

  std::unordered_map<std::string, std::size_t> state_parameter_indices = { { "A_B_eq", 0 } };
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  // State with 1 grid cell
  DenseMatrixPolicy state(1, num_species);
  state[0] = { 0.01, 0.02, 0.05 };

  // State parameters: K_eq for each grid cell (1 cell, 1 parameter)
  DenseMatrixPolicy state_parameters(1, 1, 3.3e-2);

  set.SubtractJacobianTerms(state, state_parameters, jacobian);

  // G = K_eq * [A] * [B] - [AB]
  // dG/d[A] = K_eq * [B] = 3.3e-2 * 0.02 = 0.00066
  // dG/d[B] = K_eq * [A] = 3.3e-2 * 0.01 = 0.00033
  // dG/d[AB] = -1

  // Jacobian subtracts these values (matching ProcessSet convention)
  // Constraint replaces row 2 (AB's row)
  EXPECT_NEAR(jacobian[0][2][0], -0.00066, 1e-10);  // J[2, A] -= dG/dA
  EXPECT_NEAR(jacobian[0][2][1], -0.00033, 1e-10);  // J[2, B] -= dG/dB
  EXPECT_NEAR(jacobian[0][2][2], 1.0, 1e-10);       // J[2, AB] -= dG/dAB = -(-1) = 1
}

template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testEmptyConstraintSet()
{
  // Empty constraint set should be valid and do nothing
  ConstraintSetPolicy set;

  EXPECT_EQ(set.Size(), 0);

  auto non_zero_elements = set.NonZeroJacobianElements();
  EXPECT_TRUE(non_zero_elements.empty());

  // AddForcingTerms and SubtractJacobianTerms should do nothing for empty set
  DenseMatrixPolicy state(1, 3);
  state[0] = { 0.1, 0.2, 0.3 };

  DenseMatrixPolicy forcing(1, 4, 1.0);

  // Empty state_parameters for empty constraint set
  DenseMatrixPolicy state_parameters(1, 0);

  set.AddForcingTerms(state, state_parameters, forcing);

  // Forcing should be unchanged
  EXPECT_DOUBLE_EQ(forcing[0][0], 1.0);
  EXPECT_DOUBLE_EQ(forcing[0][1], 1.0);
  EXPECT_DOUBLE_EQ(forcing[0][2], 1.0);
  EXPECT_DOUBLE_EQ(forcing[0][3], 1.0);
}

template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testUnknownSpeciesThrows()
{
  // Creating a constraint with unknown species should throw
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "invalid",
      std::vector<StoichSpecies>{ StoichSpecies(Species("X"), 1.0), StoichSpecies(Species("Y"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("XY"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 } };

  EXPECT_THROW((ConstraintSetPolicy(std::move(constraints), variable_map)), micm::MicmException);
}

/// @brief Test 3D state (3 species) with 1 constraint
template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testThreeDStateOneConstraint()
{
  const double K_eq = 3.3e-2;
  const std::size_t num_species = 3;

  // Create constraint: X <-> Y with K_eq = 3.3e-2
  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "X_Y_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("X"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("Y"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "X", 0 }, { "Y", 1 }, { "Z", 2 } };

  ConstraintSetPolicy set(std::move(constraints), variable_map);

  EXPECT_EQ(set.Size(), 1);

  // Check non-zero Jacobian elements
  // Algebraic species = Y (index 1), constraint replaces row 1
  auto non_zero_elements = set.NonZeroJacobianElements();
  EXPECT_EQ(non_zero_elements.size(), 2);                       // (1,0) dG/dX, (1,1) dG/dY + diagonal
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(1, 0)));   // dG/dX at row 1, col 0
  EXPECT_TRUE(non_zero_elements.count(std::make_pair(1, 1)));   // dG/dY at row 1, col 1 + diagonal
  EXPECT_FALSE(non_zero_elements.count(std::make_pair(1, 2)));  // Z not involved

  // Build sparse Jacobian (3x3)
  auto builder = SparseMatrixPolicy::Create(num_species)
                     .SetNumberOfBlocks(2)  // Test with 2 grid cells
                     .InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrixPolicy jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  std::unordered_map<std::string, std::size_t> state_parameter_indices = { { "X_Y_eq", 0 } };
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  // State with 2 grid cells
  DenseMatrixPolicy state(2, num_species);
  // Grid cell 0: Away from equilibrium
  state[0][0] = 10.0;  // X
  state[0][1] = 0.2;   // Y
  state[0][2] = 0.5;   // Z (uninvolved)
  // Grid cell 1: At equilibrium (Y/X = 3.3e-2)
  state[1][0] = 10.0;  // X
  state[1][1] = 0.33;  // Y = 3.3e-2 * 10.0 = 0.33
  state[1][2] = 0.3;   // Z

  // Test forcing terms
  DenseMatrixPolicy forcing(2, num_species, 0.0);

  // State parameters: K_eq for each grid cell (2 cells, 1 parameter)
  DenseMatrixPolicy state_parameters(2, 1, 3.3e-2);

  set.AddForcingTerms(state, state_parameters, forcing);

  // Constraint replaces row 1 (Y's row)
  // Grid cell 0: G = K_eq * [X] - [Y] = 3.3e-2 * 10.0 - 0.2 = 0.33 - 0.2 = 0.13
  EXPECT_NEAR(forcing[0][1], 0.13, 1e-10);
  // Grid cell 1: G = 3.3e-2 * 10.0 - 0.33 = 0.0 (at equilibrium)
  EXPECT_NEAR(forcing[1][1], 0.0, 1e-10);

  // Test Jacobian terms
  set.SubtractJacobianTerms(state, state_parameters, jacobian);

  // For constraint G = K_eq * [X] - [Y]:
  // dG/dX = K_eq = 3.3e-2
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
template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testFourDStateTwoConstraints()
{
  const double K_eq1 = 3.3e-2;
  const double K_eq2 = 3.3e-2;
  const std::size_t num_species = 4;

  // Create two constraints
  std::vector<Constraint> constraints;

  // Constraint 1: A <-> B with K_eq1 = 3.3e-2, algebraic species = B (row 1)
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  // Constraint 2: C + D <-> A with K_eq2 = 3.3e-2, algebraic species = A (row 0)
  constraints.push_back(EquilibriumConstraint(
      "CD_A_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("C"), 1.0), StoichSpecies(Species("D"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "C", 2 }, { "D", 3 } };

  ConstraintSetPolicy set(std::move(constraints), variable_map);

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
  auto builder = SparseMatrixPolicy::Create(num_species)
                     .SetNumberOfBlocks(3)  // Test with 3 grid cells
                     .InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrixPolicy jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  std::unordered_map<std::string, std::size_t> state_parameter_indices = { { "A_B_eq", 0 }, { "CD_A_eq", 1 } };
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  // State with 3 grid cells
  DenseMatrixPolicy state(3, num_species);

  // Grid cell 0: Both constraints satisfied
  state[0][0] = 0.33;     // A
  state[0][1] = 0.01089;  // B = K_eq1 * A = 3.3e-2 * 0.33
  state[0][2] = 10.0;     // C
  state[0][3] = 1.0;      // D, so C*D = 10.0, K_eq2*C*D = 3.3e-2 * 10 = 0.33 = A

  // Grid cell 1: First constraint satisfied, second not
  state[1][0] = 0.33;     // A
  state[1][1] = 0.01089;  // B = K_eq1 * A (constraint 1 satisfied)
  state[1][2] = 5.0;      // C
  state[1][3] = 1.0;      // D, C*D = 5.0, K_eq2*C*D = 3.3e-2 * 5 = 0.165 ≠ 0.33

  // Grid cell 2: Neither constraint satisfied
  state[2][0] = 0.33;  // A
  state[2][1] = 0.02;  // B ≠ K_eq1 * A = 0.01089
  state[2][2] = 5.0;   // C
  state[2][3] = 1.0;   // D, C*D = 5.0, K_eq2*C*D = 3.3e-2 * 5 = 0.165 ≠ 0.33

  // Test forcing terms
  DenseMatrixPolicy forcing(3, num_species, 0.0);

  // State parameters: K_eq for each constraint (3 cells, 2 parameters)
  DenseMatrixPolicy state_parameters(3, 2);
  for (std::size_t i = 0; i < 3; ++i)
  {
    state_parameters[i][0] = 3.3e-2;  // K_eq1 for A_B_eq
    state_parameters[i][1] = 3.3e-2;  // K_eq2 for CD_A_eq
  }

  set.AddForcingTerms(state, state_parameters, forcing);

  // Constraint 1 replaces row 1, Constraint 2 replaces row 0
  // Grid cell 0: Both at equilibrium
  EXPECT_NEAR(forcing[0][1], 0.0, 1e-5);   // G1 = K_eq1 * 0.33 - 0.01089 ≈ 0
  EXPECT_NEAR(forcing[0][0], 0.0, 1e-10);  // G2 = K_eq2 * 10.0 * 1.0 - 0.33 = 0

  // Grid cell 1: First satisfied, second not
  EXPECT_NEAR(forcing[1][1], 0.0, 1e-5);      // G1 = K_eq1 * 0.33 - 0.01089 ≈ 0
  EXPECT_NEAR(forcing[1][0], -0.165, 1e-10);  // G2 = K_eq2 * 5.0 * 1.0 - 0.33 = 0.165 - 0.33 = -0.165

  // Grid cell 2: Neither satisfied
  EXPECT_NEAR(forcing[2][1], -0.00911, 1e-5);  // G1 = K_eq1 * 0.33 - 0.02 = 0.01089 - 0.02 = -0.00911
  EXPECT_NEAR(forcing[2][0], -0.165, 1e-10);   // G2 = K_eq2 * 5.0 * 1.0 - 0.33 = -0.165

  // Test Jacobian terms
  set.SubtractJacobianTerms(state, state_parameters, jacobian);

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
template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testCoupledConstraintsSharedSpecies()
{
  const double K_eq1 = 3.3e-2;
  const double K_eq2 = 3.3e-2;
  const std::size_t num_species = 3;

  std::vector<Constraint> constraints;

  // Both constraints depend on species A
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("B"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  constraints.push_back(EquilibriumConstraint(
      "A_C_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("C"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "C", 2 } };

  ConstraintSetPolicy set(std::move(constraints), variable_map);

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
  auto builder = SparseMatrixPolicy::Create(num_species).SetNumberOfBlocks(1).InitialValue(0.0);

  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrixPolicy jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  std::unordered_map<std::string, std::size_t> state_parameter_indices = { { "A_B_eq", 0 }, { "A_C_eq", 1 } };
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  // State at dual equilibrium: [B]/[A] = 3.3e-2, [C]/[A] = 3.3e-2
  DenseMatrixPolicy state(1, num_species);
  state[0][0] = 0.1;     // A
  state[0][1] = 0.0033;  // B = K_eq1 * A = 3.3e-2 * 0.1
  state[0][2] = 0.0033;  // C = K_eq2 * A = 3.3e-2 * 0.1

  // Test forcing terms
  DenseMatrixPolicy forcing(1, num_species, 0.0);

  // State parameters: K_eq for each constraint (1 cell, 2 parameters)
  DenseMatrixPolicy state_parameters(1, 2);
  state_parameters[0][0] = 3.3e-2;  // K_eq1 for A_B_eq
  state_parameters[0][1] = 3.3e-2;  // K_eq2 for A_C_eq

  set.AddForcingTerms(state, state_parameters, forcing);

  // Both constraints should be satisfied
  EXPECT_NEAR(forcing[0][1], 0.0, 1e-10);  // G1 at row 1
  EXPECT_NEAR(forcing[0][2], 0.0, 1e-10);  // G2 at row 2

  // Test Jacobian terms
  set.SubtractJacobianTerms(state, state_parameters, jacobian);

  // Constraint 1 at row 1: dG1/dA = K_eq1 = 3.3e-2, dG1/dB = -1
  EXPECT_NEAR(jacobian[0][1][0], -K_eq1, 1e-10);
  EXPECT_NEAR(jacobian[0][1][1], 1.0, 1e-10);

  // Constraint 2 at row 2: dG2/dA = K_eq2 = 3.3e-2, dG2/dC = -1
  EXPECT_NEAR(jacobian[0][2][0], -K_eq2, 1e-10);
  EXPECT_NEAR(jacobian[0][2][2], 1.0, 1e-10);
}

template<class DenseMatrixPolicy, class SparseMatrixPolicy, class ConstraintSetPolicy>
void testVectorizedMatricesRespectGridCellIndexing()
{
  const std::size_t num_species = 3;

  std::vector<Constraint> constraints;
  constraints.push_back(EquilibriumConstraint(
      "A_B_eq",
      std::vector<StoichSpecies>{ StoichSpecies(Species("A"), 1.0), StoichSpecies(Species("B"), 1.0) },
      std::vector<StoichSpecies>{ StoichSpecies(Species("AB"), 1.0) },
      VantHoffParam{ .K_HLC_ref = 3.3e-2, .delta_H = -24000.0 }));

  std::unordered_map<std::string, std::size_t> variable_map = { { "A", 0 }, { "B", 1 }, { "AB", 2 } };

  ConstraintSetPolicy set(std::move(constraints), variable_map);
  auto non_zero_elements = set.NonZeroJacobianElements();

  // Constraint replaces AB's row (index 2), Jacobian is 3x3
  auto builder = SparseMatrixPolicy::Create(num_species).SetNumberOfBlocks(3).InitialValue(0.0);
  for (std::size_t i = 0; i < num_species; ++i)
    builder = builder.WithElement(i, i);
  for (const auto& elem : non_zero_elements)
    builder = builder.WithElement(elem.first, elem.second);

  SparseMatrixPolicy jacobian{ builder };
  set.SetJacobianFlatIds(jacobian);

  std::unordered_map<std::string, std::size_t> state_parameter_indices = { { "A_B_eq", 0 } };
  set.SetConstraintFunctions(variable_map, state_parameter_indices, jacobian);

  DenseMatrixPolicy state(3, num_species, 0.0);
  state[0] = { 0.01, 0.02, 0.05 };
  state[1] = { 0.03, 0.01, 0.2 };
  state[2] = { 0.001, 0.002, 0.004 };

  DenseMatrixPolicy forcing(3, num_species, 0.0);

  // State parameters: K_eq for each grid cell (3 cells, 1 parameter)
  DenseMatrixPolicy state_parameters(3, 1, 3.3e-2);

  set.AddForcingTerms(state, state_parameters, forcing);

  // Constraint residual replaces row 2 (AB)
  // K_eq = 3.3e-2
  EXPECT_NEAR(forcing[0][2], 3.3e-2 * 0.01 * 0.02 - 0.05, 1e-9);
  EXPECT_NEAR(forcing[1][2], 3.3e-2 * 0.03 * 0.01 - 0.2, 1e-9);
  EXPECT_NEAR(forcing[2][2], 3.3e-2 * 0.001 * 0.002 - 0.004, 1e-9);

  set.SubtractJacobianTerms(state, state_parameters, jacobian);

  // Jacobian entries at row 2 (AB's row, replaced by constraint)
  EXPECT_NEAR(jacobian[0][2][0], -(3.3e-2 * 0.02), 1e-12);
  EXPECT_NEAR(jacobian[0][2][1], -(3.3e-2 * 0.01), 1e-12);
  EXPECT_NEAR(jacobian[0][2][2], 1.0, 1e-12);

  EXPECT_NEAR(jacobian[1][2][0], -(3.3e-2 * 0.01), 1e-12);
  EXPECT_NEAR(jacobian[1][2][1], -(3.3e-2 * 0.03), 1e-12);
  EXPECT_NEAR(jacobian[1][2][2], 1.0, 1e-12);

  EXPECT_NEAR(jacobian[2][2][0], -(3.3e-2 * 0.002), 1e-12);
  EXPECT_NEAR(jacobian[2][2][1], -(3.3e-2 * 0.001), 1e-12);
  EXPECT_NEAR(jacobian[2][2][2], 1.0, 1e-12);
}
