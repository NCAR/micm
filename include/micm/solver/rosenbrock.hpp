/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 * 
 * Much of this solver was formulated and implemented from this book:
 * Hairer, E., Wanner, G., 1996. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, 2nd edition. ed. Springer, Berlin ; New York.
 * The source code for many (all?) of the solvers in that book can be found here: http://www.unige.ch/~hairer/software.html
 * 
 * Some extensions to the rosenbrock solver formulated there were formulated in this paper
 * Sandu, A., Verwer, J.G., Blom, J.G., Spee, E.J., Carmichael, G.R., Potra, F.A., 1997. Benchmarking stiff ode solvers for atmospheric chemistry problems II: Rosenbrock solvers. Atmospheric Environment 31, 3459–3472. https://doi.org/10.1016/S1352-2310(97)83212-8
 * 
 */
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <iostream>
#include <limits>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/solver.hpp>
#include <string>
#include <vector>

namespace micm
{

  /**
   * @brief An implementation of the Chapman mechnanism solver
   *
   */
  class RosenbrockSolver : public Solver
  {
   private:
    struct Rosenbrock_params
    {
      size_t N_{};
      size_t stages_{};
      size_t upper_limit_tolerance_{};
      size_t max_number_of_steps_{100};

      double round_off_{std::numeric_limits<double>::epsilon()}; // Unit roundoff (1+round_off)>1
      double factor_min_{0.2};                 // solver step size minimum boundary
      double factor_max_{6};                 // solver step size maximum boundary
      double rejection_factor_decrease_{0.1};  // used to decrease the step after 2 successive rejections
      double safety_factor_{0.9};              // safety factor in new step size computation

      double h_min_{0};    // step size min
      double h_max_{0.5};    // step size max
      double h_start_{0.005};  // step size start

      std::array<bool, 6>
          new_function_evaluation_{};  // which steps reuse the previous iterations evaluation or do a new evaluation

      double estimator_of_local_order_{};  // the minumum between the main and the embedded scheme orders plus one
      std::array<double, 15> a_{};         // coefficient matrix a
      std::array<double, 15> c_{};         // coefficient matrix c
      std::array<double, 6> m_{};          // coefficients for new step evaluation
      std::array<double, 6> e_{};          // error estimation coefficients
      std::array<double, 6> alpha_{};
      std::array<double, 6> gamma_{};

      double absolute_tolerance_{ 1e-12 };
      double relative_tolerance_{ 1e-4 };
    };

   public:
    Rosenbrock_params parameters_;
    std::vector<double> rate_constants_;
    Solver::Rosenbrock_stats stats_;
    static constexpr uint64_t number_sparse_factor_elements_ = 23;

    static constexpr double delta_min_ = 1.0e-5;
    static constexpr double error_min_ = 1.0e-10;

   public:
    /// @brief Default constructor
    RosenbrockSolver();
    ~RosenbrockSolver();

    /// @brief A virtual function to be defined by any solver baseclass
    /// @param time_start Time step to start at
    /// @param time_end Time step to end at
    /// @param number_densities Species concentrations in molecules / cm3
    /// @param number_density_air The number density of air in molecules / cm3
    /// @return A struct containing results and a status code
    SolverResult Solve(
        const double& time_start,
        const double& time_end,
        const std::vector<double>& number_densities,
        const double& number_density_air) noexcept override;

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
    std::vector<double> force(
        const std::vector<double>& rate_constants,
        const std::vector<double>& number_densities,
        const double& number_density_air);

    /// @brief compute jacobian decomposition of [alpha * I - dforce_dy]
    /// @param dforce_dy
    /// @param alpha
    /// @return An jacobian decomposition
    std::vector<double> factored_alpha_minus_jac(const std::vector<double>& dforce_dy, const double& alpha);

    /// @brief Computes product of [dforce_dy * vector]
    /// @param dforce_dy  jacobian of forcing
    /// @param vector vector ordered as the order of number density in dy
    /// @return Product of jacobian with vector
    std::vector<double> dforce_dy_times_vector(const std::vector<double>& dforce_dy, const std::vector<double>& vector);

    /// @brief Update the rate constants for the environment state
    /// @param temperature in kelvin
    /// @param pressure in pascals
    void calculate_rate_constants(const double& temperature, const double& pressure);

    /// @brief Solve the system
    /// @param K idk, something
    /// @param ode_jacobian the jacobian
    /// @return the new state?
    std::vector<double> lin_solve(const std::vector<double>& K, const std::vector<double>& ode_jacobian);

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
        const double& number_density_air);

   private:
    /// @brief Factor
    /// @param jacobian
    void factor(std::vector<double>& jacobian);

    std::vector<double> backsolve_L_y_eq_b(const std::vector<double>& jacobian, const std::vector<double>& b);
    std::vector<double> backsolve_U_x_eq_b(const std::vector<double>& jacobian, const std::vector<double>& y);

    /// @brief Initializes the solving parameters for a three-stage rosenbrock solver
    void three_stage_rosenbrock();


    /// @brief Computes the scaled norm of the vector errors
    /// @param original_number_densities the original number densities
    /// @param new_number_densities the new number densities
    /// @param errors The computed errors
    /// @return
    double error_norm(
        std::vector<double> original_number_densities,
        std::vector<double> new_number_densities,
        std::vector<double> errors);
  };

  inline RosenbrockSolver::RosenbrockSolver()
      : parameters_(),
        rate_constants_(7, 0),
        stats_()
  {
    three_stage_rosenbrock();
  }

  inline RosenbrockSolver::~RosenbrockSolver()
  {
  }

  inline Solver::SolverResult RosenbrockSolver::Solve(
      const double& time_start,
      const double& time_end,
      const std::vector<double>& original_number_densities,
      const double& number_density_air) noexcept
  {
    std::vector<std::vector<double>> K(parameters_.stages_, std::vector<double>(parameters_.N_, 0));
    std::vector<double> Y(original_number_densities);
    std::vector<double> Ynew(original_number_densities.size(), 0);
    std::vector<double> forcing{};

    double present_time = time_start;
    double H =
        std::min(std::max(std::abs(parameters_.h_min_), std::abs(parameters_.h_start_)), std::abs(parameters_.h_max_));

    SolverResult result{};
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

      auto initial_forcing = force(rate_constants_, Y, number_density_air);

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
        auto ode_jacobian = lin_factor(H, parameters_.gamma_[0], is_singular, Y, number_density_air);
        stats_.jacobian_updates += 1;
        if (is_singular)
        {
          result.state_ = Solver::SolverState::RepeatedlySingularMatrix;
          break;
        }

        // Compute the stages
        {
          // the first stage (stage 0), inlined to remove a branch in the following for loop
          K[0] = lin_solve(initial_forcing, ode_jacobian);

          // stages (1-# of stages)
          for (uint64_t stage = 1; stage < parameters_.stages_; ++stage)
          {
            double stage_combinations = ((stage+1) - 1) * ((stage+1) - 2) / 2;
            if (parameters_.new_function_evaluation_[stage])
            {
              Ynew = Y;
              for (uint64_t j = 0; j <= stage; ++j)
              {
                auto a = parameters_.a_[stage_combinations + j];
                for (uint64_t idx = 0; idx < Ynew.size(); ++idx)
                {
                  Ynew[idx] += a * K[j][idx];
                }
              }
              forcing = force(rate_constants_, Ynew, number_density_air);
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
            K[stage] = lin_solve(K[stage], ode_jacobian);
          }
        }

        // Compute the new solution
        Ynew = Y;
        for (uint64_t stage = 0; stage < parameters_.stages_; ++stage)
        {
          for (uint64_t idx = 0; idx < Ynew.size(); ++idx)
          {
            Ynew[idx] += parameters_.m_[stage] * K[stage][idx];
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
        auto error = error_norm(Y, Ynew, Yerror);

        // New step size is bounded by FacMin <= Hnew/H <= FacMax
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
          Y = Ynew;
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

  inline std::vector<std::string> RosenbrockSolver::reaction_names()
  {
  }

  inline std::vector<std::string> RosenbrockSolver::photolysis_names()
  {
  }

  inline std::vector<std::string> RosenbrockSolver::species_names()
  {
  }

  inline std::vector<double> RosenbrockSolver::force(
      const std::vector<double>& rate_constants,
      const std::vector<double>& number_densities,
      const double& number_density_air)
  {
    std::vector<double> force(number_densities.size(), 0);

    stats_.function_calls += 1;
    return force;
  }

  inline std::vector<double> RosenbrockSolver::factored_alpha_minus_jac(
      const std::vector<double>& dforce_dy,
      const double& alpha)
  {
    std::vector<double> jacobian(number_sparse_factor_elements_);

    factor(jacobian);
    return jacobian;
  }

  inline void RosenbrockSolver::factor(std::vector<double>& jacobian)
  {
  }

  inline std::vector<double> RosenbrockSolver::dforce_dy_times_vector(
      const std::vector<double>& dforce_dy,
      const std::vector<double>& vector)
  {
    std::vector<double> result(dforce_dy.size(), 0);

    return result;
  }

  inline std::vector<double> RosenbrockSolver::backsolve_L_y_eq_b(
      const std::vector<double>& jacobian,
      const std::vector<double>& b)
  {
    std::vector<double> y(parameters_.N_, 0);
    return y;
  }

  inline std::vector<double> RosenbrockSolver::backsolve_U_x_eq_b(
      const std::vector<double>& jacobian,
      const std::vector<double>& y)
  {
    std::vector<double> x(y.size(), 0);
    double temporary{};
    return x;
  }

  inline void RosenbrockSolver::three_stage_rosenbrock()
  {
    // an L-stable method, 3 stages, order 3, 2 function evaluations
    //
    // original formaulation for three stages:
    // Sandu, A., Verwer, J.G., Blom, J.G., Spee, E.J., Carmichael, G.R., Potra, F.A., 1997. 
    // Benchmarking stiff ode solvers for atmospheric chemistry problems II: Rosenbrock solvers. 
    // Atmospheric Environment 31, 3459–3472. https://doi.org/10.1016/S1352-2310(97)83212-8


    parameters_.stages_ = 3;
    parameters_.N_ = species_names().size();

    //  The coefficient matrices A and C are strictly lower triangular.
    //  The lower triangular (subdiagonal) elements are stored in row-wise order:
    //  A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
    //  The general mapping formula is:
    //      A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
    //      C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    parameters_.a_.fill(0);
    parameters_.a_[0] = 1;
    parameters_.a_[1] = 1;
    parameters_.a_[2] = 0;

    parameters_.c_.fill(0);
    parameters_.c_[0] = -0.10156171083877702091975600115545e+01;
    parameters_.c_[1] = 0.40759956452537699824805835358067e+01;
    parameters_.c_[2] = 0.92076794298330791242156818474003e+01;

    // Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
    // or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    parameters_.new_function_evaluation_.fill(false);
    parameters_.new_function_evaluation_[0] = true;
    parameters_.new_function_evaluation_[1] = true;
    parameters_.new_function_evaluation_[2] = false;

    // Coefficients for new step solution
    parameters_.m_.fill(0);
    parameters_.m_[0] = 0.1e+01;
    parameters_.m_[1] = 0.61697947043828245592553615689730e+01;
    parameters_.m_[2] = -0.42772256543218573326238373806514;

    // Coefficients for error estimator
    parameters_.e_.fill(0);
    parameters_.e_[0] = 0.5;
    parameters_.e_[1] = -0.29079558716805469821718236208017e+01;
    parameters_.e_[2] = 0.22354069897811569627360909276199;

    // ros_ELO = estimator of local order - the minimum between the
    // main and the embedded scheme orders plus 1
    parameters_.estimator_of_local_order_ = 3;

    // Y_stage_i ~ Y( T + H*Alpha_i )
    parameters_.alpha_.fill(0);
    parameters_.alpha_[0] = 0;
    parameters_.alpha_[1] = 0.43586652150845899941601945119356;
    parameters_.alpha_[2] = 0.43586652150845899941601945119356;

    // Gamma_i = \sum_j  gamma_{i,j}
    parameters_.gamma_.fill(0);
    parameters_.gamma_[0] = 0.43586652150845899941601945119356;
    parameters_.gamma_[1] = 0.24291996454816804366592249683314;
    parameters_.gamma_[2] = 0.21851380027664058511513169485832e+01;
  }

  inline void RosenbrockSolver::calculate_rate_constants(const double& temperature, const double& pressure)
  {
  }

  inline std::vector<double> RosenbrockSolver::lin_factor(
      double& H,
      const double& gamma,
      bool& singular,
      const std::vector<double>& number_densities,
      const double& number_density_air)
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
      ode_jacobian = factored_alpha_minus_jac(dforce_dy(rate_constants_, number_densities, number_density_air), alpha);
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

  inline std::vector<double> RosenbrockSolver::lin_solve(const std::vector<double>& K, const std::vector<double>& jacobian)
  {
    auto y = backsolve_L_y_eq_b(jacobian, K);
    auto x = backsolve_U_x_eq_b(jacobian, y);
    stats_.solves += 1;
    return x;
  }

  inline std::vector<double> RosenbrockSolver::dforce_dy(
      const std::vector<double>& rate_constants,
      const std::vector<double>& number_densities,
      const double& number_density_air)
  {
    std::vector<double> jacobian(number_sparse_factor_elements_, 0);

    return jacobian;
  }

  inline double RosenbrockSolver::error_norm(
      std::vector<double> Y,
      std::vector<double> Ynew,
      std::vector<double> errors)
  {
    // Solving Ordinary Differential Equations II, page 123
    // https://link-springer-com.cuucar.idm.oclc.org/book/10.1007/978-3-642-05221-7
    std::vector<double> maxs(Y.size());
    std::vector<double> scale(Y.size());

    for (uint64_t idx = 0; idx < Y.size(); ++idx)
    {
      maxs[idx] = std::max(std::abs(Y[idx]), std::abs(Ynew[idx]));
    }

    for (uint64_t idx = 0; idx < Y.size(); ++idx)
    {
      scale[idx] = parameters_.absolute_tolerance_ + parameters_.relative_tolerance_ * maxs[idx];
    }

    double sum = 0;
    for (uint64_t idx = 0; idx < Y.size(); ++idx)
    {
      sum += std::pow(errors[idx] / scale[idx], 2);
    }

    return std::max(std::sqrt(sum / parameters_.N_), error_min_);
  }
}  // namespace micm