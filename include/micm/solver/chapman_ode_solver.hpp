/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

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
    void Solve(double state[]) override;

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
    std::vector<double> p_force(std::vector<double> rate_constants, std::vector<double> number_densities, double number_density_air);
  };

  inline ChapmanODESolver::ChapmanODESolver()
  {
  }

  inline ChapmanODESolver::~ChapmanODESolver()
  {
  }

  inline void ChapmanODESolver::Solve(double state[])
  {
  }

  inline std::vector<std::string> ChapmanODESolver::reaction_names()
  {
    return std::vector<std::string>{ "O2_1", "O3_1", "O3_2", "N2_O1D_1", "O1D_O2_1", "O_O3_1", "M_O_O2_1" };
  }

  inline std::vector<std::string> ChapmanODESolver::photolysis_names()
  {
    return std::vector<std::string>{ "O2_1", "O3_1", "O3_2", };
  }

  inline std::vector<std::string> ChapmanODESolver::species_names()
  {
    return std::vector<std::string>{ "M", "Ar", "CO2", "H2O", "N2", "O1D", "O", "O2", "O3", };
  }

  inline std::vector<double> ChapmanODESolver::p_force(std::vector<double> rate_constants, std::vector<double> number_densities, double number_density_air)
  {
    return std::vector<double>();
  }

}  // namespace micm