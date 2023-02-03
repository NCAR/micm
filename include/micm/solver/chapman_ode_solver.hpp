/* Copyright (C) 2023 National Center for Atmospheric Research,
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <micm/solver/solver.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <string>
#include <vector>
#include <array>

namespace micm
{

  /**
   * @brief An implementation of the Chapman mechnanism solver
   *
   */
  class ChapmanODESolver : public Solver
  {
  private:
    struct Rosenbrock_params {
      size_t N_ {};
      size_t stages_ {};
      size_t upper_limit_tolerance_ {};
      size_t max_number_of_steps_ {};

      double round_off_ {}; // Unit roundoff (1+round_off)>1
      double factor_min_ {}; // solver step size minimum boundary
      double factor_max_ {}; // solver step size maximum boundary
      double rejection_factor_decrease_ {}; // used to decrease the step after 2 successive rejections
      double safety_factor_ {}; // safety factor in new step size computation

      double h_min_ {}; // step size min
      double h_max_ {}; // step size max
      double h_start_ {}; // step size start

      std::array<bool, 6> new_function_evaluation_ {}; // which steps reuse the previous iterations evaluation or do a new evaluation

      double estimator_of_local_order_ {}; // the minumu between the main and the embedded scheme orders plus one
      std::array<double, 15> a_ {}; // coefficient matrix a
      std::array<double, 15> c_ {}; // coefficient matrix c
      std::array<double, 6> m_ {}; // coefficients for new step evaluation
      std::array<double, 6> e_ {}; // error estimation coefficients
      std::array<double, 6> alpha_ {};
      std::array<double, 6> gamma_ {};

      std::vector<double> absolute_tolerance_ {};
      std::vector<double> relative_tolerance_ {};
    };

    struct Rosenbrock_stats {
      uint64_t forcing_function_calls {}; // Nfun
      uint64_t jacobian_updates {}; // Njac
      uint64_t number_of_steps {}; // Nstp
      uint64_t accepted {}; // Nacc
      uint64_t rejected {}; // Nrej
      uint64_t decompositions {}; // Ndec
      uint64_t solves {}; // Nsol
      uint64_t singular {}; // Nsng
      uint64_t total_steps {}; // Ntotstp

      void reset(){
        forcing_function_calls = 0;
        jacobian_updates = 0;
        number_of_steps = 0;
        accepted = 0;
        rejected = 0;
        decompositions = 0;
        solves = 0;
        singular = 0;
        total_steps = 0;
      }

    };

  public:
    Rosenbrock_params parameters_;
    Rosenbrock_stats stats_;
    std::array<double, 7> rate_constants_;

    static constexpr double delta_min_ = 1.0e-5;

   public:
    /// @brief Default constructor
    ChapmanODESolver();
    ~ChapmanODESolver();

    /// @brief An implementation of the 3-stage rosenbrock solver
    /// @param time_start Time step to start at
    /// @param time_end Time step to end at
    /// @param number_densities Species concentrations in molecules / cm3
    /// @return A SolverResult containing the new state and solver status codes
    SolverResult Solve(double time_start, double time_end, std::vector<double> number_densities) override;

    std::vector<double> matrix_solver(std::vector<double> LU, std::vector<double> b);

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

    /// @brief Update the rate constants for the environment state
    /// @param temperature in kelvin
    /// @param pressure in pascals
    void calculate_rate_constants(double temperature, double pressure);

   private:
    /// @brief Factor
    /// @param LU
    void factor(std::vector<double>& LU);

    std::vector<double> backsolve_L_y_eq_b(std::vector<double>& LU, std::vector<double>& b);
    std::vector<double> backsolve_U_x_eq_b(std::vector<double>& LU, std::vector<double>& y);

    /// @brief Initializes the solving parameters for a three-stage rosenbrock solver
    void three_stage_rosenbrock();
  };

  inline ChapmanODESolver::ChapmanODESolver()
    : parameters_(), 
      stats_(),
      rate_constants_()
  {
    three_stage_rosenbrock();
  }

  inline ChapmanODESolver::~ChapmanODESolver()
  {
  }

  inline Solver::SolverResult ChapmanODESolver::Solve(double time_start, double time_end, std::vector<double> number_densities)
  {
    double present_time = time_start;
    double H = std::min(std::max(std::abs(parameters_.h_min_),abs(parameters_.h_start_)) , abs(parameters_.h_max_) );

    SolverResult result {};

    if(std::abs(H) <= 10*parameters_.round_off_) {
      H = delta_min_;
    }

    bool reject_last_h = false;
    bool reject_more_h = false;


    // TimeLoop: 
    // DO WHILE ( (presentTime-time_end)+this%Roundoff <= ZERO )
    while( (present_time - time_end + parameters_.round_off_) <= 0){
      if (stats_.number_of_steps > parameters_.max_number_of_steps_){
        // too many steps
        result.state_ = Solver::SolverState::ConvergenceExceededMaxSteps;
        break;
      }

      if ( ((present_time + 0.1*H) == present_time) || (H <= parameters_.round_off_) ) {
        result.state_ = Solver::SolverState::StepSizeTooSmall;
        break;
      }

      //  Limit H if necessary to avoid going beyond time_end
      H = std::min(H,std::abs(time_end-present_time));

      //   Compute the function at current time
      // Fcn0(:) = p_force( Y )
      // this%icntrl(Nfun) = this%icntrl(Nfun) + 1


    // ======= left off here =============

    //   !~~~>  Repeat step calculation until current step accepted
    //   UntilAccepted: 
    //   DO
    //     !~~~>  Form and factor the rosenbrock ode jacobian
    //     CALL theKinetics%LinFactor( H, this%ros_Gamma(1), Y, Singular, this%icntrl )
    //     this%icntrl(Njac) = this%icntrl(Njac) + 1
    //     IF (Singular) THEN ! More than 5 consecutive failed decompositions
    //         Ierr = -8
    //         CALL ros_ErrorMsg(-8,present_time,H,IERR)
    //         RETURN
    //     END IF

    //     !~~~>   Compute the stages
    //     Stage_loop: &
    //     DO istage = 1, this%ros_S
    //       IF ( istage /= 1 ) THEN
    //         S_ndx = (istage - 1)*(istage - 2)/2
    //         IF ( this%ros_NewF(istage) ) THEN
    //           Ynew(1:N) = Y(1:N)
    //           DO j = 1, istage-1
    //             Ynew(1:N) = Ynew(1:N) + this%ros_A(S_ndx+j)*K(1:N,j)
    //           END DO
    //           Tau = present_time + this%ros_Alpha(istage)*H
    //           Fcn(:) = theKinetics%force( Ynew )
    //           this%icntrl(Nfun) = this%icntrl(Nfun) + 1
    //         ENDIF
    //         K(:,istage) = Fcn(:)
    //         DO j = 1, istage-1
    //           HC = this%ros_C(S_ndx+j)/H
    //           K(1:N,istage) = K(1:N,istage) + HC*K(1:N,j)
    //         END DO
    //       ELSE
    //         K(:,1) = Fcn0(:)
    //         Fcn(:) = Fcn0(:)
    //       ENDIF
    //       CALL theKinetics%LinSolve( K(:,istage) )
    //       this%icntrl(Nsol) = this%icntrl(Nsol) + 1
    //     END DO Stage_loop

    //     !~~~>  Compute the new solution
    //     Ynew(1:N) = Y(1:N)
    //     DO j=1,this%ros_S
    //       Ynew(1:N) = Ynew(1:N) + this%ros_M(j)*K(1:N,j)
    //     END DO

    //     !~~~>  Compute the error estimation
    //     Yerr(1:N) = ZERO
    //     DO j=1,this%ros_S
    //       Yerr(1:N) = Yerr(1:N) + this%ros_E(j)*K(1:N,j)
    //     END DO
    //     Err = ros_ErrorNorm( this, Y, Ynew, Yerr )

    //     !~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
    //     Fac  = MIN(this%FacMax,MAX(this%FacMin,this%FacSafe/Err**(ONE/this%ros_ELO)))
    //     Hnew = H*Fac

    //     !~~~>  Check the error magnitude and adjust step size
    //     this%icntrl(Nstp) = this%icntrl(Nstp) + 1
    //     this%icntrl(Ntotstp) = this%icntrl(Ntotstp) + 1
    //     Accepted: &
    //     IF ( (Err <= ONE).OR.(H <= this%Hmin) ) THEN
    //       this%icntrl(Nacc) = this%icntrl(Nacc) + 1
    //       Y(1:N) = Ynew(1:N)
    //       present_time = present_time + H
    //       Hnew = MAX(this%Hmin,MIN(Hnew,this%Hmax))
    //       IF (RejectLastH) THEN  ! No step size increase after a rejected step
    //           Hnew = MIN(Hnew,H)
    //       END IF
    //       this%rcntrl(Nhexit) = H
    //       this%rcntrl(Nhnew)  = Hnew
    //       this%rcntrl(Ntexit) = present_time
    //       RejectLastH = .FALSE.
    //       RejectMoreH = .FALSE.
    //       H = Hnew
    //       EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
    //     ELSE Accepted  !~~~> Reject step
    //       IF (RejectMoreH) THEN
    //           Hnew = H*this%FacRej
    //       END IF
    //       RejectMoreH = RejectLastH
    //       RejectLastH = .TRUE.
    //       H = Hnew
    //       IF (this%icntrl(Nacc) >= 1)  this%icntrl(Nrej) = this%icntrl(Nrej) + 1
    //     ENDIF Accepted ! Err <= 1

    //   END DO UntilAccepted
    }

    // !~~~> Succesful exit
    // IERR = 0  !~~~> The integration was successful
    // IF( present(T) ) THEN
    //   T = presentTime
    // ENDIF

    return result;
  }

  inline std::vector<double> ChapmanODESolver::matrix_solver(std::vector<double> LU, std::vector<double> b)
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

  inline void ChapmanODESolver::three_stage_rosenbrock()
  {
   // an L-stable method, 3 stages, order 3, 2 function evaluations

    parameters_.stages_ = 3;

   //  The coefficient matrices A and C are strictly lower triangular.
   //  The lower triangular (subdiagonal) elements are stored in row-wise order:
   //  A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
   //  The general mapping formula is:
   //      A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
   //      C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   parameters_.a_.fill(0);
   parameters_.a_[0] = 1;
   parameters_.a_[0] = 1;
   parameters_.a_[0] = 0;

   parameters_.c_.fill(0);
   parameters_.c_[0] = -0.10156171083877702091975600115545e+01;
   parameters_.c_[0] = 0.40759956452537699824805835358067e+01;
   parameters_.c_[0] = 0.92076794298330791242156818474003e+01;

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

  inline void ChapmanODESolver::calculate_rate_constants(double temperature, double pressure){
    //O2_1
    //k_O2_1: O2 -> 2*O
    // photolysis = rate_constant_photolysis_t( &
    //     photolysis_rate_constant_index = 1 )
    // rate_constants(1) = photolysis%calculate( environment )
    rate_constants_[0] = 0; // TODO fix

    //O3_1
    //k_O3_1: O3 -> 1*O1D + 1*O2
    // photolysis = rate_constant_photolysis_t( &
    //     photolysis_rate_constant_index = 2 )
    // rate_constants(2) = photolysis%calculate( environment )
    rate_constants_[1] = 0; // TODO fix

    //O3_2
    //k_O3_2: O3 -> 1*O + 1*O2
    // photolysis = rate_constant_photolysis_t( &
    //     photolysis_rate_constant_index = 3 )
    // rate_constants(3) = photolysis%calculate( environment )
    rate_constants_[2] = 0; // TODO fix

    //N2_O1D_1
    //k_N2_O1D_1: N2 + O1D -> 1*O + 1*N2
    rate_constants_[3] = ArrheniusRateConstant(2.15e-11, 0, 110, 0, 0).calculate(temperature, pressure);

    //O1D_O2_1
    //k_O1D_O2_1: O1D + O2 -> 1*O + 1*O2
    rate_constants_[4] = ArrheniusRateConstant(3.3e-11, 0, 55, 0, 0).calculate(temperature, pressure);

    //O_O3_1
    //k_O_O3_1: O + O3 -> 2*O2
    rate_constants_[5] = ArrheniusRateConstant(8e-12, 0, -2060, 0, 0).calculate(temperature, pressure);

    //M_O_O2_1
    //k_M_O_O2_1: M + O + O2 -> 1*O3 + 1*M
    rate_constants_[5] = ArrheniusRateConstant(6e-34, 0, 2.4, 0, 0).calculate(temperature, pressure);
  }

}  // namespace micm