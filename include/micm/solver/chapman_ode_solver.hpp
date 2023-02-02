/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <micm/solver/solver.hpp>
#include <string>
#include <vector>

namespace micm
{

  /**
   * @brief An implementation of the Chapman mechnanism solver
   *
   */
  class ChapmanODESolver : public Solver
  {
   private:
   public:
    /// @brief Default constructor
    ChapmanODESolver();
    ~ChapmanODESolver();
    /// @brief Move the system to the next state
    /// @param state The collection of species concentrations
    std::vector<double> Solve(std::vector<double> LU, std::vector<double> b) override;

    /// @brief Returns a list of reaction names
    /// @return vector of strings
    std::vector<std::string> reaction_names();

    /// @brief Returns a list of species that participate in photolysis
    /// @return vector of strings
    std::vector<std::string> photolysis_names();

    /// @brief Returns a list of species names
    /// @return vector of strings
    std::vector<std::string> species_names();

    /// @brief Calculate a chemical forcing
    /// @param rate_constants List of rate constants for each needed species
    /// @param number_densities The number density of each species
    /// @param number_density_air The number density of air
    /// @return A vector of forcings
    std::vector<double>
    p_force(std::vector<double> rate_constants, std::vector<double> number_densities, double number_density_air);

    /// @brief compute LU decomposition of [alpha * I - dforce_dy]
    /// @param dforce_dy
    /// @param alpha
    /// @return An LU decomposition
    std::vector<double> factored_alpha_minus_jac(std::vector<double> dforce_dy, double alpha);

    /// @brief Computes product of [dforce_dy * vector]
    /// @param dforce_dy  Jacobian of forcing
    /// @param vector vector ordered as the order of number density in dy
    /// @return Product of jacobian with vector
    std::vector<double> dforce_dy_times_vector(std::vector<double> dforce_dy, std::vector<double> vector);

   private:
    /// @brief Factor
    /// @param LU
    void factor(std::vector<double>& LU);

    std::vector<double> backsolve_L_y_eq_b(std::vector<double>& LU, std::vector<double>& b);
    std::vector<double> backsolve_U_x_eq_b(std::vector<double>& LU, std::vector<double>& y);
  };

  inline ChapmanODESolver::ChapmanODESolver()
  {
  }

  inline ChapmanODESolver::~ChapmanODESolver()
  {
  }

  inline std::vector<double> ChapmanODESolver::Solve(std::vector<double> LU, std::vector<double> b)
  {
    auto y = backsolve_L_y_eq_b(LU, b);
    auto x = backsolve_U_x_eq_b(LU, y);
    return x;
  }

  inline std::vector<std::string> ChapmanODESolver::reaction_names()
  {
    return std::vector<std::string>{ "O2_1", "O3_1", "O3_2", "N2_O1D_1", "O1D_O2_1", "O_O3_1", "M_O_O2_1" };
  }

  inline std::vector<std::string> ChapmanODESolver::photolysis_names()
  {
    return std::vector<std::string>{
      "O2_1",
      "O3_1",
      "O3_2",
    };
  }

  inline std::vector<std::string> ChapmanODESolver::species_names()
  {
    return std::vector<std::string>{
      "M", "Ar", "CO2", "H2O", "N2", "O1D", "O", "O2", "O3",
    };
  }

  inline std::vector<double> ChapmanODESolver::p_force(
      std::vector<double> rate_constants,
      std::vector<double> number_densities,
      double number_density_air)
  {
    // Forcings:
    // M, Ar, CO2, H2O, N2, O1D, O, O2, O3,
    std::vector<double> force(number_densities.size(), 0);

    assert(force.size() == 9);

    // M, Ar, CO2, H2O, N2 are all zero

    // O1D
    {
      // k_O3_1: O3 -> 1*O1D + 1*O2
      force[5] = force[5] + rate_constants[1] * number_densities[8];
      // k_N2_O1D_1: N2 + O1D -> 1*O + 1*N2
      force[5] = force[5] - rate_constants[3] * number_densities[4] * number_densities[5];
      // k_O1D_O2_1: O1D + O2 -> 1*O + 1*O2
      force[5] = force[5] - rate_constants[4] * number_densities[5] * number_densities[7];
    }

    // O
    {
      // k_O2_1: O2 -> 2*O
      force[6] = force[6] + 2 * rate_constants[0] * number_densities[7];
      // k_O3_2: O3 -> 1*O + 1*O2
      force[6] = force[6] + rate_constants[2] * number_densities[8];
      // k_N2_O1D_1: N2 + O1D -> 1*O + 1*N2
      force[6] = force[6] + rate_constants[3] * number_densities[4] * number_densities[5];
      // k_O1D_O2_1: O1D + O2 -> 1*O + 1*O2
      force[6] = force[6] + rate_constants[4] * number_densities[5] * number_densities[7];
      // k_O_O3_1: O + O3 -> 2*O2
      force[6] = force[6] - rate_constants[5] * number_densities[6] * number_densities[8];
      // k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
      force[6] = force[6] - rate_constants[6] * number_densities[0] * number_densities[6] * number_densities[7];
    }

    // O2
    {
      // k_O2_1: O2 -> 2*O
      force[7] = force[7] - rate_constants[0] * number_densities[7];
      // k_O3_1: O3 -> 1*O1D + 1*O2
      force[7] = force[7] + rate_constants[1] * number_densities[8];
      // k_O3_2: O3 -> 1*O + 1*O2
      force[7] = force[7] + rate_constants[2] * number_densities[8];
      // k_O_O3_1: O + O3 -> 2*O2
      force[7] = force[7] + 2 * rate_constants[5] * number_densities[6] * number_densities[8];
      // k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
      force[7] = force[7] - rate_constants[6] * number_densities[0] * number_densities[6] * number_densities[7];
    }

    // O3
    {
      // k_O3_1: O3 -> 1*O1D + 1*O2
      force[8] = force[8] - rate_constants[1] * number_densities[8];
      // k_O3_2: O3 -> 1*O + 1*O2
      force[8] = force[8] - rate_constants[2] * number_densities[8];
      // k_O_O3_1: O + O3 -> 2*O2
      force[8] = force[8] - rate_constants[5] * number_densities[6] * number_densities[8];
      // k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
      force[8] = force[8] + rate_constants[6] * number_densities[0] * number_densities[6] * number_densities[7];
    }

    return force;
  }

  inline std::vector<double> ChapmanODESolver::factored_alpha_minus_jac(std::vector<double> dforce_dy, double alpha)
  {
    std::vector<double> LU(dforce_dy);
    // multiply LU by -1
    std::transform(LU.begin(), LU.end(), LU.begin(), [](auto& c) { return -c; });

    assert(LU.size() >= 23);

    LU[0] = -dforce_dy[0] + alpha;
    LU[4] = -dforce_dy[4] + alpha;
    LU[5] = -dforce_dy[5] + alpha;
    LU[6] = -dforce_dy[6] + alpha;
    LU[7] = -dforce_dy[7] + alpha;
    LU[10] = -dforce_dy[10] + alpha;
    LU[12] = -dforce_dy[12] + alpha;
    LU[17] = -dforce_dy[17] + alpha;
    LU[22] = -dforce_dy[22] + alpha;

    factor(LU);
    return LU;
  }

  inline void ChapmanODESolver::factor(std::vector<double>& LU)
  {
    LU[0] = 1. / LU[0];
    LU[1] = LU[1] * LU[0];
    LU[2] = LU[2] * LU[0];
    LU[3] = LU[3] * LU[0];
    LU[4] = 1. / LU[4];
    LU[5] = 1. / LU[5];
    LU[6] = 1. / LU[6];
    LU[7] = 1. / LU[7];
    LU[8] = LU[8] * LU[7];
    LU[9] = LU[9] * LU[7];
    LU[10] = 1. / LU[10];
    LU[11] = LU[11] * LU[10];
    LU[16] = LU[16] - LU[11] * LU[15];
    LU[20] = LU[20] - LU[11] * LU[19];
    LU[12] = 1. / LU[12];
    LU[13] = LU[13] * LU[12];
    LU[14] = LU[14] * LU[12];
    LU[17] = LU[17] - LU[13] * LU[16];
    LU[18] = LU[18] - LU[14] * LU[16];
    LU[21] = LU[21] - LU[13] * LU[20];
    LU[22] = LU[22] - LU[14] * LU[20];
    LU[17] = 1. / LU[17];
    LU[18] = LU[18] * LU[17];
    LU[22] = LU[22] - LU[18] * LU[21];
    LU[22] = 1. / LU[22];
  }

  inline std::vector<double> ChapmanODESolver::dforce_dy_times_vector(
      std::vector<double> dforce_dy,
      std::vector<double> vector)
  {
    std::vector<double> result(dforce_dy.size(), 0);

    assert(result.size() >= 23);

    // df_O/d[M] * M_temporary
    result[6] = result[6] + dforce_dy[1] * vector[0];
    // df_O2/d[M] * M_temporary
    result[7] = result[7] + dforce_dy[2] * vector[0];
    // df_O3/d[M] * M_temporary
    result[8] = result[8] + dforce_dy[3] * vector[0];
    // df_O1D/d[N2] * N2_temporary
    result[5] = result[5] + dforce_dy[8] * vector[4];
    // df_O/d[N2] * N2_temporary
    result[6] = result[6] + dforce_dy[9] * vector[4];
    // df_O1D/d[O1D] * O1D_temporary
    result[5] = result[5] + dforce_dy[10] * vector[5];
    // df_O/d[O1D] * O1D_temporary
    result[6] = result[6] + dforce_dy[11] * vector[5];
    // df_O/d[O] * O_temporary
    result[6] = result[6] + dforce_dy[12] * vector[6];
    // df_O2/d[O] * O_temporary
    result[7] = result[7] + dforce_dy[13] * vector[6];
    // df_O3/d[O] * O_temporary
    result[8] = result[8] + dforce_dy[14] * vector[6];
    // df_O1D/d[O2] * O2_temporary
    result[5] = result[5] + dforce_dy[15] * vector[7];
    // df_O/d[O2] * O2_temporary
    result[6] = result[6] + dforce_dy[16] * vector[7];
    // df_O2/d[O2] * O2_temporary
    result[7] = result[7] + dforce_dy[17] * vector[7];
    // df_O3/d[O2] * O2_temporary
    result[8] = result[8] + dforce_dy[18] * vector[7];
    // df_O1D/d[O3] * O3_temporary
    result[5] = result[5] + dforce_dy[19] * vector[8];
    // df_O/d[O3] * O3_temporary
    result[6] = result[6] + dforce_dy[20] * vector[8];
    // df_O2/d[O3] * O3_temporary
    result[7] = result[7] + dforce_dy[21] * vector[8];
    // df_O3/d[O3] * O3_temporary
    result[8] = result[8] + dforce_dy[22] * vector[8];

    return result;
  }


  inline std::vector<double> ChapmanODESolver::backsolve_L_y_eq_b(std::vector<double>& LU, std::vector<double>& b){
    std::vector<double> y(LU.size());

    y[0] = b[0];
    y[1] = b[1];
    y[2] = b[2];
    y[3] = b[3];
    y[4] = b[4];
    y[5] = b[5];
    y[5] = y[5] - LU[8] * y[4];
    y[6] = b[6];
    y[6] = y[6] - LU[1] * y[0];
    y[6] = y[6] - LU[9] * y[4];
    y[6] = y[6] - LU[11] * y[5];
    y[7] = b[7];
    y[7] = y[7] - LU[2] * y[0];
    y[7] = y[7] - LU[13] * y[6];
    y[8] = b[8];
    y[8] = y[8] - LU[3] * y[0];
    y[8] = y[8] - LU[14] * y[6];
    y[8] = y[8] - LU[18] * y[7];

    return y;
  }

  inline std::vector<double> ChapmanODESolver::backsolve_U_x_eq_b(std::vector<double>& LU, std::vector<double>& y){
    std::vector<double> x(LU.size(), 0);
    double temporary{};

    temporary = y[8];
    x[8] = LU[22] * temporary;
    temporary = y[7];
    temporary = temporary - LU[21] * x[8];
    x[7] = LU[17] * temporary;
    temporary = y[6];
    temporary = temporary - LU[16] * x[7];
    temporary = temporary - LU[20] * x[8];
    x[6] = LU[12] * temporary;
    temporary = y[5];
    temporary = temporary - LU[15] * x[7];
    temporary = temporary - LU[19] * x[8];
    x[5] = LU[10] * temporary;
    temporary = y[4];
    x[4] = LU[7] * temporary;
    temporary = y[3];
    x[3] = LU[6] * temporary;
    temporary = y[2];
    x[2] = LU[5] * temporary;
    temporary = y[1];
    x[1] = LU[4] * temporary;
    temporary = y[0];
    x[0] = LU[0] * temporary;

    return x;
  }

}  // namespace micm