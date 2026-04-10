// Copyright (C) 2024-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>

namespace micm
{

  /// Result of comparing an analytical Jacobian against a finite-difference approximation
  struct JacobianComparisonResult
  {
    bool passed{ true };
    double max_abs_error{ 0.0 };
    double max_rel_error{ 0.0 };
    std::size_t worst_block{ 0 };
    std::size_t worst_row{ 0 };
    std::size_t worst_col{ 0 };
    double worst_analytical{ 0.0 };
    double worst_fd{ 0.0 };
  };

  /// Compute a dense finite-difference Jacobian approximation using central differences.
  ///
  /// The forcing callable should have the signature:
  ///   void(const DenseMatrixPolicy& variables, DenseMatrixPolicy& forcing)
  /// where variables has shape [num_blocks x num_species] and forcing has the same shape.
  /// Callers should bind any additional arguments (rate constants, state parameters, etc.)
  /// into the callable via a lambda capture.
  ///
  /// Returns a DenseMatrixPolicy of shape [num_blocks x (num_species * num_species)] where
  /// element [block][row * num_species + col] = df_row/dx_col.
  ///
  /// When a perturbation would push a variable below zero, one-sided differences are used instead.
  template<class DenseMatrixPolicy, class ForcingFunc>
  DenseMatrixPolicy FiniteDifferenceJacobian(
      ForcingFunc forcing_func,
      const DenseMatrixPolicy& base_variables,
      std::size_t num_species,
      double perturbation = 1.0e-8)
  {
    const std::size_t num_blocks = base_variables.NumRows();
    DenseMatrixPolicy result(num_blocks, num_species * num_species, 0.0);

    // Workspace matrices
    DenseMatrixPolicy vars_plus(num_blocks, num_species, 0.0);
    DenseMatrixPolicy vars_minus(num_blocks, num_species, 0.0);
    DenseMatrixPolicy forcing_plus(num_blocks, num_species, 0.0);
    DenseMatrixPolicy forcing_minus(num_blocks, num_species, 0.0);
    DenseMatrixPolicy forcing_base(num_blocks, num_species, 0.0);

    for (std::size_t col = 0; col < num_species; ++col)
    {
      // Set up perturbed states
      vars_plus.Copy(base_variables);
      vars_minus.Copy(base_variables);

      for (std::size_t block = 0; block < num_blocks; ++block)
      {
        double x_j = base_variables[block][col];
        double h = perturbation * std::max(1.0, std::abs(x_j));
        vars_plus[block][col] = x_j + h;
        vars_minus[block][col] = x_j - h;
      }

      forcing_plus.Fill(0.0);
      forcing_minus.Fill(0.0);

      forcing_func(vars_plus, forcing_plus);
      forcing_func(vars_minus, forcing_minus);

      for (std::size_t block = 0; block < num_blocks; ++block)
      {
        double x_j = base_variables[block][col];
        double h = perturbation * std::max(1.0, std::abs(x_j));

        // Check if the negative perturbation would go below zero
        bool use_central = (x_j - h >= 0.0);

        if (use_central)
        {
          double two_h = 2.0 * h;
          for (std::size_t row = 0; row < num_species; ++row)
          {
            result[block][row * num_species + col] = (forcing_plus[block][row] - forcing_minus[block][row]) / two_h;
          }
        }
        else
        {
          // Forward difference: need f(x) as well
          forcing_base.Fill(0.0);
          forcing_func(base_variables, forcing_base);
          for (std::size_t row = 0; row < num_species; ++row)
          {
            result[block][row * num_species + col] = (forcing_plus[block][row] - forcing_base[block][row]) / h;
          }
        }
      }
    }

    return result;
  }

  /// Compare an analytical sparse Jacobian (which stores -df/dx per MICM convention)
  /// against a finite-difference dense Jacobian (which stores +df/dx).
  ///
  /// Uses a combined tolerance: |analytical - fd| < atol + rtol * max(|analytical|, |fd|)
  ///
  /// Only non-zero elements in the sparse Jacobian are compared.
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  JacobianComparisonResult CompareJacobianToFiniteDifference(
      const SparseMatrixPolicy& analytical_jacobian,
      const DenseMatrixPolicy& fd_jacobian,
      std::size_t num_species,
      double atol = 1.0e-7,
      double rtol = 1.0e-7)
  {
    JacobianComparisonResult result;
    const std::size_t num_blocks = analytical_jacobian.NumberOfBlocks();

    for (std::size_t block = 0; block < num_blocks; ++block)
    {
      for (std::size_t row = 0; row < num_species; ++row)
      {
        for (std::size_t col = 0; col < num_species; ++col)
        {
          if (analytical_jacobian.IsZero(row, col))
            continue;

          // Analytical stores -(df/dx); negate to get +df/dx
          double analytical_val = -analytical_jacobian[block][row][col];
          double fd_val = fd_jacobian[block][row * num_species + col];

          double abs_error = std::abs(analytical_val - fd_val);
          double scale = std::max(std::abs(analytical_val), std::abs(fd_val));
          double rel_error = (scale > 0.0) ? abs_error / scale : 0.0;

          if (abs_error > result.max_abs_error)
          {
            result.max_abs_error = abs_error;
            result.worst_block = block;
            result.worst_row = row;
            result.worst_col = col;
            result.worst_analytical = analytical_val;
            result.worst_fd = fd_val;
          }
          if (rel_error > result.max_rel_error)
          {
            result.max_rel_error = rel_error;
          }

          double tolerance = atol + rtol * scale;
          if (abs_error > tolerance)
          {
            result.passed = false;
          }
        }
      }
    }

    return result;
  }

  /// Check that no significant Jacobian entry exists outside the declared sparsity pattern.
  ///
  /// This catches missing NonZeroJacobianElements declarations. Any FD entry outside the
  /// sparsity pattern that exceeds the threshold indicates an undeclared dependency.
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  JacobianComparisonResult CheckJacobianSparsityCompleteness(
      const SparseMatrixPolicy& analytical_jacobian,
      const DenseMatrixPolicy& fd_jacobian,
      std::size_t num_species,
      double threshold = 1.0e-6)
  {
    JacobianComparisonResult result;
    const std::size_t num_blocks = analytical_jacobian.NumberOfBlocks();

    for (std::size_t block = 0; block < num_blocks; ++block)
    {
      for (std::size_t row = 0; row < num_species; ++row)
      {
        for (std::size_t col = 0; col < num_species; ++col)
        {
          if (!analytical_jacobian.IsZero(row, col))
            continue;

          double fd_val = std::abs(fd_jacobian[block][row * num_species + col]);
          if (fd_val > threshold)
          {
            result.passed = false;
            if (fd_val > result.max_abs_error)
            {
              result.max_abs_error = fd_val;
              result.worst_block = block;
              result.worst_row = row;
              result.worst_col = col;
              result.worst_analytical = 0.0;
              result.worst_fd = fd_val;
            }
          }
        }
      }
    }

    return result;
  }

}  // namespace micm
