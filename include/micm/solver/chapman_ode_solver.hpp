/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <iostream>
#include <limits>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <string>
#include <vector>

namespace micm
{

  /**
   * @brief An implementation of the Chapman mechnanism solver
   *
   */
  class ChapmanODESolver : public RosenbrockSolver<micm::Matrix>
  {
    std::size_t number_sparse_factor_elements_ = 23;
   public:
    /// @brief Default constructor
    ChapmanODESolver();
    ~ChapmanODESolver();

    /// Returns a state variable for the Chapman system
    /// @return State variable for Chapman
    State<> GetState() const;

    /// @brief A virtual function to be defined by any solver baseclass
    /// @param time_start Time step to start at
    /// @param time_end Time step to end at
    /// @param state The system state to solve for
    /// @return A struct containing results and a status code
    Solver::SolverResult Solve(
        double time_start,
        double time_end,
        State<>& state) noexcept;

    /// @brief Returns a list of reaction names
    /// @return vector of strings
    std::vector<std::string> reaction_names() override;

    /// @brief Returns a list of species that participate in photolysis
    /// @return vector of strings
    std::vector<std::string> photolysis_names() override;

    /// @brief Returns a list of species names
    /// @return vector of strings
    std::vector<std::string> species_names() override;

    /// @brief Calculate a chemical forcing
    /// @param rate_constants List of rate constants for each needed species
    /// @param number_densities The number density of each species
    /// @param number_density_air The number density of air
    /// @return A vector of forcings
    std::vector<double> force(
        const std::vector<double>& rate_constants,
        const std::vector<double>& number_densities,
        const double& number_density_air);

    /// @brief compute jacobian decomposition of [alpha * I - dforce_dy]
    /// @param dforce_dy
    /// @param alpha
    /// @return An jacobian decomposition
    std::vector<double> factored_alpha_minus_jac(const std::vector<double>& dforce_dy, const double& alpha) override;

    /// @brief Computes product of [dforce_dy * vector]
    /// @param dforce_dy  jacobian of forcing
    /// @param vector vector ordered as the order of number density in dy
    /// @return Product of jacobian with vector
    std::vector<double> dforce_dy_times_vector(const std::vector<double>& dforce_dy, const std::vector<double>& vector)
        override;

    /// @brief Update the rate constants for the environment state
    /// @param state The current state of the chemical system
    void UpdateState(State<>& state);

    /// @brief Solve the system
    /// @param K idk, something
    /// @param ode_jacobian the jacobian
    /// @return the new state?
    std::vector<double> lin_solve(const std::vector<double>& K, const std::vector<double>& ode_jacobian) override;

    /// @brief Compute the derivative of the forcing w.r.t. each chemical, the jacobian
    /// @param rate_constants List of rate constants for each needed species
    /// @param number_densities The number density of each species
    /// @param number_density_air The number density of air
    /// @return The jacobian
    std::vector<double> dforce_dy(
        const std::vector<double>& rate_constants,
        const std::vector<double>& number_densities,
        const double& number_density_air);

    /// @brief Prepare the rosenbrock ode solver matrix
    /// @param H time step (seconds)
    /// @param gamma time step factor for specific rosenbrock method
    /// @param Y  constituent concentration (molec/cm^3)
    /// @param singular indicates if the matrix is singular
    std::vector<double> lin_factor(
        double& H,
        const double& gamma,
        bool& singular,
        const std::vector<double>& number_densities,
        const double& number_density_air,
        const std::vector<double>& rate_constants);

    /// @brief Factor
    /// @param jacobian
    void factor(std::vector<double>& jacobian) override;

    std::vector<double> backsolve_L_y_eq_b(const std::vector<double>& jacobian, const std::vector<double>& b) override;
    std::vector<double> backsolve_U_x_eq_b(const std::vector<double>& jacobian, const std::vector<double>& y) override;
  };

  inline ChapmanODESolver::ChapmanODESolver()
      : RosenbrockSolver()
  {
    three_stage_rosenbrock();
    // override state size for hard-coded Chapman mechanism
    parameters_.N_ = 9;
  }

  inline ChapmanODESolver::~ChapmanODESolver()
  {
  }

  inline State<> ChapmanODESolver::GetState() const
  {
    return State{ 9, 3, 7 };
  }

  inline Solver::SolverResult ChapmanODESolver::Solve(
      double time_start,
      double time_end,
      State<>& state) noexcept
  {
    std::vector<std::vector<double>> K(parameters_.stages_, std::vector<double>(parameters_.N_, 0));
    std::vector<double> Y(state.variables_[0]);
    std::vector<double> rate_constants = state.rate_constants_[0];
    const double number_density_air = state.conditions_[0].air_density_;
    std::vector<double> forcing{};

    double present_time = time_start;
    double H =
        std::min(std::max(std::abs(parameters_.h_min_), std::abs(parameters_.h_start_)), std::abs(parameters_.h_max_));

    Solver::SolverResult result{};
    stats_.reset();

    if (std::abs(H) <= 10 * parameters_.round_off_)
    {
      H = delta_min_;
    }

    bool reject_last_h = false;
    bool reject_more_h = false;

    while ((present_time - time_end + parameters_.round_off_) <= 0)
    {
      if (stats_.number_of_steps > parameters_.max_number_of_steps_)
      {
        result.state_ = Solver::SolverState::ConvergenceExceededMaxSteps;
        break;
      }

      if (((present_time + 0.1 * H) == present_time) || (H <= parameters_.round_off_))
      {
        result.state_ = Solver::SolverState::StepSizeTooSmall;
        break;
      }

      //  Limit H if necessary to avoid going beyond time_end
      H = std::min(H, std::abs(time_end - present_time));

      auto initial_forcing = force(rate_constants, Y, number_density_air);

      bool accepted = false;
      //  Repeat step calculation until current step accepted
      while (!accepted)
      {
        if (stats_.number_of_steps > parameters_.max_number_of_steps_)
        {
          break;
        }
        bool is_singular{ false };
        // Form and factor the rosenbrock ode jacobian
        auto ode_jacobian = lin_factor(H, parameters_.gamma_[0], is_singular, Y, number_density_air, rate_constants);
        stats_.jacobian_updates += 1;
        if (is_singular)
        {
          result.state_ = Solver::SolverState::RepeatedlySingularMatrix;
          break;
        }

        // Compute the stages
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
        {
          if (stage == 0)
          {
            K[0] = initial_forcing;
            forcing = initial_forcing;
          }
          else
          {
            double stage_combinations = ((stage + 1) - 1) * ((stage + 1) - 2) / 2;
            if (parameters_.new_function_evaluation_[stage])
            {
              auto new_Y(Y);
              for (uint64_t j = 0; j <= stage; ++j)
              {
                assert((stage_combinations + j) < parameters_.a_.size());
                auto a = parameters_.a_[stage_combinations + j];
                for (uint64_t idx = 0; idx < new_Y.size(); ++idx)
                {
                  new_Y[idx] += a * K[j][idx];
                }
              }
              forcing = force(rate_constants, new_Y, number_density_air);
            }
            K[stage] = forcing;
            for (uint64_t j = 0; j < stage; ++j)
            {
              auto HC = parameters_.c_[stage_combinations + j] / H;
              for (uint64_t idx = 0; idx < K[stage].size(); ++idx)
              {
                K[stage][idx] += HC * K[j][idx];
              }
            }
          }
          K[stage] = lin_solve(K[stage], ode_jacobian);
        }

        // Compute the new solution
        auto new_Y(Y);
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
        {
          for (uint64_t idx = 0; idx < new_Y.size(); ++idx)
          {
            new_Y[idx] += parameters_.m_[stage] * K[stage][idx];
          }
        }

        // Compute the error estimation
        std::vector<double> Yerror(Y.size(), 0);
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
        {
          uint64_t idx = 0;
          for (uint64_t idx = 0; idx < Yerror.size(); ++idx)
          {
            Yerror[idx] += parameters_.e_[stage] * K[stage][idx];
          }
        }
        auto error = error_norm(Y, new_Y, Yerror);

        // New step size is bounded by FacMin <= Hnew/H <= FacMax
        // Fac  = MIN(this%FacMax,MAX(this%FacMin,this%FacSafe/Err**(ONE/this%ros_ELO)))
        double Hnew = H * std::min(
                              parameters_.factor_max_,
                              std::max(
                                  parameters_.factor_min_,
                                  parameters_.safety_factor_ / std::pow(error, 1 / parameters_.estimator_of_local_order_)));

        // Check the error magnitude and adjust step size
        stats_.number_of_steps += 1;
        stats_.total_steps += 1;
        if ((error < 1) || (H < parameters_.h_min_))
        {
          stats_.accepted += 1;
          present_time = present_time + H;
          Y = new_Y;
          Hnew = std::max(parameters_.h_min_, std::min(Hnew, parameters_.h_max_));
          if (reject_last_h)
          {
            // No step size increase after a rejected step
            Hnew = std::min(Hnew, H);
          }
          reject_last_h = false;
          reject_more_h = false;
          H = Hnew;
          accepted = true;
        }
        else
        {
          // Reject step
          if (reject_more_h)
          {
            Hnew = H * parameters_.rejection_factor_decrease_;
          }
          reject_more_h = reject_last_h;
          reject_last_h = true;
          H = Hnew;
          if (stats_.accepted >= 1)
          {
            stats_.rejected += 1;
          }
        }
      }
    }

    result.T = present_time;
    result.stats_ = stats_;
    result.result_ = std::move(Y);
    result.state_ = Solver::SolverState::Converged;

    return result;
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

  inline std::vector<double> ChapmanODESolver::force(
      const std::vector<double>& rate_constants,
      const std::vector<double>& number_densities,
      const double& number_density_air)
  {
    // Forcings:
    // M, Ar, CO2, H2O, N2, O1D, O, O2, O3,
    std::vector<double> force(number_densities.size(), 0);

    assert(force.size() >= 9);

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

    stats_.function_calls += 1;
    return force;
  }

  inline std::vector<double> ChapmanODESolver::factored_alpha_minus_jac(
      const std::vector<double>& dforce_dy,
      const double& alpha)
  {
    std::vector<double> jacobian(number_sparse_factor_elements_);
    // multiply jacobian by -1
    std::transform(dforce_dy.begin(), dforce_dy.end(), jacobian.begin(), [](auto& c) { return -c; });

    assert(jacobian.size() >= 23);

    jacobian[0] = -dforce_dy[0] + alpha;
    jacobian[4] = -dforce_dy[4] + alpha;
    jacobian[5] = -dforce_dy[5] + alpha;
    jacobian[6] = -dforce_dy[6] + alpha;
    jacobian[7] = -dforce_dy[7] + alpha;
    jacobian[10] = -dforce_dy[10] + alpha;
    jacobian[12] = -dforce_dy[12] + alpha;
    jacobian[17] = -dforce_dy[17] + alpha;
    jacobian[22] = -dforce_dy[22] + alpha;

    factor(jacobian);
    return jacobian;
  }

  inline void ChapmanODESolver::factor(std::vector<double>& jacobian)
  {
    jacobian[0] = 1. / jacobian[0];
    jacobian[1] = jacobian[1] * jacobian[0];
    jacobian[2] = jacobian[2] * jacobian[0];
    jacobian[3] = jacobian[3] * jacobian[0];
    jacobian[4] = 1. / jacobian[4];
    jacobian[5] = 1. / jacobian[5];
    jacobian[6] = 1. / jacobian[6];
    jacobian[7] = 1. / jacobian[7];
    jacobian[8] = jacobian[8] * jacobian[7];
    jacobian[9] = jacobian[9] * jacobian[7];
    jacobian[10] = 1. / jacobian[10];
    jacobian[11] = jacobian[11] * jacobian[10];
    jacobian[16] = jacobian[16] - jacobian[11] * jacobian[15];
    jacobian[20] = jacobian[20] - jacobian[11] * jacobian[19];
    jacobian[12] = 1. / jacobian[12];
    jacobian[13] = jacobian[13] * jacobian[12];
    jacobian[14] = jacobian[14] * jacobian[12];
    jacobian[17] = jacobian[17] - jacobian[13] * jacobian[16];
    jacobian[18] = jacobian[18] - jacobian[14] * jacobian[16];
    jacobian[21] = jacobian[21] - jacobian[13] * jacobian[20];
    jacobian[22] = jacobian[22] - jacobian[14] * jacobian[20];
    jacobian[17] = 1. / jacobian[17];
    jacobian[18] = jacobian[18] * jacobian[17];
    jacobian[22] = jacobian[22] - jacobian[18] * jacobian[21];
    jacobian[22] = 1. / jacobian[22];
  }

  inline std::vector<double> ChapmanODESolver::dforce_dy_times_vector(
      const std::vector<double>& dforce_dy,
      const std::vector<double>& vector)
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

  inline std::vector<double> ChapmanODESolver::backsolve_L_y_eq_b(
      const std::vector<double>& jacobian,
      const std::vector<double>& b)
  {
    std::vector<double> y(parameters_.N_, 0);

    y[0] = b[0];
    y[1] = b[1];
    y[2] = b[2];
    y[3] = b[3];
    y[4] = b[4];
    y[5] = b[5];
    y[5] = y[5] - jacobian[8] * y[4];
    y[6] = b[6];
    y[6] = y[6] - jacobian[1] * y[0];
    y[6] = y[6] - jacobian[9] * y[4];
    y[6] = y[6] - jacobian[11] * y[5];
    y[7] = b[7];
    y[7] = y[7] - jacobian[2] * y[0];
    y[7] = y[7] - jacobian[13] * y[6];
    y[8] = b[8];
    y[8] = y[8] - jacobian[3] * y[0];
    y[8] = y[8] - jacobian[14] * y[6];
    y[8] = y[8] - jacobian[18] * y[7];

    return y;
  }

  inline std::vector<double> ChapmanODESolver::backsolve_U_x_eq_b(
      const std::vector<double>& jacobian,
      const std::vector<double>& y)
  {
    std::vector<double> x(y.size(), 0);
    double temporary{};

    temporary = y[8];
    x[8] = jacobian[22] * temporary;
    temporary = y[7];
    temporary = temporary - jacobian[21] * x[8];
    x[7] = jacobian[17] * temporary;
    temporary = y[6];
    temporary = temporary - jacobian[16] * x[7];
    temporary = temporary - jacobian[20] * x[8];
    x[6] = jacobian[12] * temporary;
    temporary = y[5];
    temporary = temporary - jacobian[15] * x[7];
    temporary = temporary - jacobian[19] * x[8];
    x[5] = jacobian[10] * temporary;
    temporary = y[4];
    x[4] = jacobian[7] * temporary;
    temporary = y[3];
    x[3] = jacobian[6] * temporary;
    temporary = y[2];
    x[2] = jacobian[5] * temporary;
    temporary = y[1];
    x[1] = jacobian[4] * temporary;
    temporary = y[0];
    x[0] = jacobian[0] * temporary;

    return x;
  }

  inline void ChapmanODESolver::UpdateState(State<>& state)
  {
    double temperature = state.conditions_[0].temperature_;
    double pressure = state.conditions_[0].pressure_;

    // O2_1
    // k_O2_1: O2 -> 2*O
    state.rate_constants_[0][0] = state.custom_rate_parameters_[0][0];

    // O3_1
    // k_O3_1: O3 -> 1*O1D + 1*O2
    state.rate_constants_[0][1] = state.custom_rate_parameters_[0][1];

    // O3_2
    // k_O3_2: O3 -> 1*O + 1*O2
    state.rate_constants_[0][2] = state.custom_rate_parameters_[0][2];

    // N2_O1D_1
    // k_N2_O1D_1: N2 + O1D -> 1*O + 1*N2
    ArrheniusRateConstantParameters params;
    params.A_ = 2.15e-11;
    params.C_ = 110;
    state.rate_constants_[0][3] = ArrheniusRateConstant(params).calculate(temperature, pressure);

    // O1D_O2_1
    // k_O1D_O2_1: O1D + O2 -> 1*O + 1*O2
    params.A_ = 3.3e-11;
    params.C_ = 55;
    state.rate_constants_[0][4] = ArrheniusRateConstant(params).calculate(temperature, pressure);

    // O_O3_1
    // k_O_O3_1: O + O3 -> 2*O2
    params.A_ = 8e-12;
    params.C_ = -2060;
    state.rate_constants_[0][5] = ArrheniusRateConstant(params).calculate(temperature, pressure);

    // M_O_O2_1
    // k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
    params.A_ = 6e-34;
    params.B_ = 2.4;
    params.C_ = 0;
    state.rate_constants_[0][6] = ArrheniusRateConstant(params).calculate(temperature, pressure);
  }

  inline std::vector<double> ChapmanODESolver::lin_factor(
      double& H,
      const double& gamma,
      bool& singular,
      const std::vector<double>& number_densities,
      const double& number_density_air,
      const std::vector<double>& rate_constants)
  {
    /*
    TODO: invesitage this function. The fortran equivalent appears to have a bug.

    From my understanding the fortran do loop would only ever do one iteration and is equivalent to what's below
    */

    std::function<bool(const std::vector<double>)> is_successful = [](const std::vector<double>& jacobian) { return true; };
    std::vector<double> ode_jacobian;
    uint64_t n_consecutive = 0;
    singular = true;

    while (true)
    {
      double alpha = 1 / (H * gamma);
      // compute jacobian decomposition of alpha*I - dforce_dy
      ode_jacobian = factored_alpha_minus_jac(dforce_dy(rate_constants, number_densities, number_density_air), alpha);
      stats_.decompositions += 1;

      if (is_successful(ode_jacobian))
      {
        singular = false;
        break;
      }
      else
      {
        stats_.singular += 1;
        n_consecutive += 1;

        if (n_consecutive <= 5)
        {
          H /= 2;
        }
        else
        {
          break;
        }
      }
    }

    return ode_jacobian;
  }

  inline std::vector<double> ChapmanODESolver::lin_solve(const std::vector<double>& K, const std::vector<double>& jacobian)
  {
    auto y = backsolve_L_y_eq_b(jacobian, K);
    auto x = backsolve_U_x_eq_b(jacobian, y);
    stats_.solves += 1;
    return x;
  }

  inline std::vector<double> ChapmanODESolver::dforce_dy(
      const std::vector<double>& rate_constants,
      const std::vector<double>& number_densities,
      const double& number_density_air)
  {
    std::vector<double> jacobian(number_sparse_factor_elements_, 0);

    // df_O/d[M]
    //  k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
    jacobian[1] = jacobian[1] - rate_constants[6] * number_densities[6] * number_densities[7];

    // df_O2/d[M]
    //  k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
    jacobian[2] = jacobian[2] - rate_constants[6] * number_densities[6] * number_densities[7];

    // df_O3/d[M]
    //  k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
    jacobian[3] = jacobian[3] + rate_constants[6] * number_densities[6] * number_densities[7];

    // df_O1D/d[N2]
    //  k_N2_O1D_1: N2 + O1D -> 1*O + 1*N2
    jacobian[8] = jacobian[8] - rate_constants[3] * number_densities[5];

    // df_O/d[N2]
    //  k_N2_O1D_1: N2 + O1D -> 1*O + 1*N2
    jacobian[9] = jacobian[9] + rate_constants[3] * number_densities[5];

    // df_O1D/d[O1D]
    //  k_N2_O1D_1: N2 + O1D -> 1*O + 1*N2
    jacobian[10] = jacobian[10] - rate_constants[3] * number_densities[4];

    //  k_O1D_O2_1: O1D + O2 -> 1*O + 1*O2
    jacobian[10] = jacobian[10] - rate_constants[4] * number_densities[7];

    // df_O/d[O1D]
    //  k_N2_O1D_1: N2 + O1D -> 1*O + 1*N2
    jacobian[11] = jacobian[11] + rate_constants[3] * number_densities[4];

    //  k_O1D_O2_1: O1D + O2 -> 1*O + 1*O2
    jacobian[11] = jacobian[11] + rate_constants[4] * number_densities[7];

    // df_O/d[O]
    //  k_O_O3_1: O + O3 -> 2*O2
    jacobian[12] = jacobian[12] - rate_constants[5] * number_densities[8];

    //  k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
    jacobian[12] = jacobian[12] - rate_constants[6] * number_densities[0] * number_densities[7];

    // df_O2/d[O]
    //  k_O_O3_1: O + O3 -> 2*O2
    jacobian[13] = jacobian[13] + 2 * rate_constants[5] * number_densities[8];

    //  k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
    jacobian[13] = jacobian[13] - rate_constants[6] * number_densities[0] * number_densities[7];

    // df_O3/d[O]
    //  k_O_O3_1: O + O3 -> 2*O2
    jacobian[14] = jacobian[14] - rate_constants[5] * number_densities[8];

    //  k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
    jacobian[14] = jacobian[14] + rate_constants[6] * number_densities[0] * number_densities[7];

    // df_O1D/d[O2]
    //  k_O1D_O2_1: O1D + O2 -> 1*O + 1*O2
    jacobian[15] = jacobian[15] - rate_constants[4] * number_densities[5];

    // df_O/d[O2]
    //  k_O2_1: O2 -> 2*O
    jacobian[16] = jacobian[16] + 2 * rate_constants[0];

    //  k_O1D_O2_1: O1D + O2 -> 1*O + 1*O2
    jacobian[16] = jacobian[16] + rate_constants[4] * number_densities[5];

    //  k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
    jacobian[16] = jacobian[16] - rate_constants[6] * number_densities[0] * number_densities[6];

    // df_O2/d[O2]
    //  k_O2_1: O2 -> 2*O
    jacobian[17] = jacobian[17] - rate_constants[0];

    //  k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
    jacobian[17] = jacobian[17] - rate_constants[6] * number_densities[0] * number_densities[6];

    // df_O3/d[O2]
    //  k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
    jacobian[18] = jacobian[18] + rate_constants[6] * number_densities[0] * number_densities[6];

    // df_O1D/d[O3]
    //  k_O3_1: O3 -> 1*O1D + 1*O2
    jacobian[19] = jacobian[19] + rate_constants[1];

    // df_O/d[O3]
    //  k_O3_2: O3 -> 1*O + 1*O2
    jacobian[20] = jacobian[20] + rate_constants[2];

    //  k_O_O3_1: O + O3 -> 2*O2
    jacobian[20] = jacobian[20] - rate_constants[5] * number_densities[6];

    // df_O2/d[O3]
    //  k_O3_1: O3 -> 1*O1D + 1*O2
    jacobian[21] = jacobian[21] + rate_constants[1];

    //  k_O3_2: O3 -> 1*O + 1*O2
    jacobian[21] = jacobian[21] + rate_constants[2];

    //  k_O_O3_1: O + O3 -> 2*O2
    jacobian[21] = jacobian[21] + 2 * rate_constants[5] * number_densities[6];

    // df_O3/d[O3]
    //  k_O3_1: O3 -> 1*O1D + 1*O2
    jacobian[22] = jacobian[22] - rate_constants[1];

    //  k_O3_2: O3 -> 1*O + 1*O2
    jacobian[22] = jacobian[22] - rate_constants[2];

    //  k_O_O3_1: O + O3 -> 2*O2
    jacobian[22] = jacobian[22] - rate_constants[5] * number_densities[6];

    return jacobian;
  }
}  // namespace micm