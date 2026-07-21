// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Exactly solvable quadratic differential equation (QDE) benchmarks after
// Bacsi & Kocsis, "Exactly solvable quadratic differential equation systems
// through generalized inversion" (arXiv:2206.03728; Qual. Theory Dyn. Syst.
// 22, 2023). Their B-transformation y = x / (x^T B x) maps a solvable class
// of quadratic systems onto linear ones, giving closed-form solutions
// (their eq. 17) for genuinely nonlinear, coupled quadratic ODEs -- the same
// problem class as at-most-bimolecular mass-action kinetics. Errors here are
// measured against the closed form directly: no reference solver in the loop.
//
//   QDE-2 -- their eq. (19), 2 variables:
//       x1' = -x1^2 + 4 x1 x2 - x1 + 2 x2
//       x2' =  x1^2 + 2 x2^2 + x1
//     B = [[2,1],[1,-4]], lambda = -1, M = [[0,2],[1,1]], w = (0, 1/2).
//
//   QDE-3 -- their eq. (23), 3 variables:
//       x1' =  x1^2 + 7 x2^2 - 4 x1 x2 + 5 x1
//       x2' = -2 x2^2 + x1 x2 + 2 x1 x3 - 7 x2 x3 + 2 x2
//       x3' = -x2^2 - 7 x3^2 - 4 x2 x3 - x3
//     B = [[0,0,1/2],[0,1,0],[1/2,0,0]], lambda = 4, M = diag(1,-2,-5),
//     w = (7, 2, -1).
//
// Both are run in shifted coordinates X = x + c with constant c chosen so the
// trajectory stays in the positive orthant (MICM clamps species nonnegative);
// the shifted system is quadratic + linear + constant with the same Jacobian,
// and X_exact(t) = x_exact(t) + c. Initial conditions keep the trajectory
// clear of the hypersurface x^T B x = 0 where the inversion degenerates (and
// where generic QDE trajectories blow up in finite time).
//
// Outputs:
//   qde_exact.csv      -- adaptive tolerance sweep vs the closed form
//   qde_fixed_step.csv -- fixed-step error at t = T for the three tableaus
#include <micm/CPU.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
  // ---- QDE-2 ----------------------------------------------------------------
  constexpr double Q2_SHIFT[2] = { 2.0, 2.0 };
  constexpr double Q2_X0[2] = { -1.2, -0.8 };  // unshifted IC
  constexpr double Q2_TEND = 3.0;

  /// Original-coordinate RHS and Jacobian of QDE-2.
  inline void Qde2Rhs(const double x[2], double f[2])
  {
    f[0] = -x[0] * x[0] + 4.0 * x[0] * x[1] - x[0] + 2.0 * x[1];
    f[1] = x[0] * x[0] + 2.0 * x[1] * x[1] + x[0];
  }

  /// Exact shifted solution of QDE-2 via Bacsi-Kocsis eq. (17).
  /// e^{Mt} for M = [[0,2],[1,1]] (eigenvalues 2, -1) in closed form.
  std::array<double, 2> Qde2Exact(double t)
  {
    const double a = std::exp(2.0 * t), b = std::exp(-t);
    const double E11 = (a + 2.0 * b) / 3.0, E12 = 2.0 * (a - b) / 3.0;
    const double E21 = (a - b) / 3.0, E22 = (2.0 * a + b) / 3.0;
    // y_p = (e^{Mt} - I) M^{-1} w with M^{-1} w = (1/2, 0)
    const double yp1 = 0.5 * (E11 - 1.0), yp2 = 0.5 * E21;
    const double B[2][2] = { { 2.0, 1.0 }, { 1.0, -4.0 } };
    const double lambda = -1.0;
    const double x01 = Q2_X0[0], x02 = Q2_X0[1];
    const double b0 = B[0][0] * x01 * x01 + 2.0 * B[0][1] * x01 * x02 + B[1][1] * x02 * x02;
    const double ex1 = E11 * x01 + E12 * x02, ex2 = E21 * x01 + E22 * x02;
    const double yBex = yp1 * (B[0][0] * ex1 + B[0][1] * ex2) + yp2 * (B[1][0] * ex1 + B[1][1] * ex2);
    const double yBy = yp1 * (B[0][0] * yp1 + B[0][1] * yp2) + yp2 * (B[1][0] * yp1 + B[1][1] * yp2);
    const double den = std::exp(-lambda * t) + 2.0 * yBex + b0 * yBy;
    return { (ex1 + b0 * yp1) / den + Q2_SHIFT[0], (ex2 + b0 * yp2) / den + Q2_SHIFT[1] };
  }

  /// QDE-2 in shifted coordinates X = x + c as a MICM external model.
  class Qde2Model
  {
   public:
    std::set<std::string> SpeciesUsed() const
    {
      return { "X1", "X2" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto x1 = s.at("X1"), x2 = s.at("X2");
      return { { x1, x1 }, { x1, x2 }, { x2, x1 }, { x2, x2 } };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto x1 = s.at("X1"), x2 = s.at("X2");
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          const double x[2] = { state[cell][x1] - Q2_SHIFT[0], state[cell][x2] - Q2_SHIFT[1] };
          double f[2];
          Qde2Rhs(x, f);
          forcing[cell][x1] += f[0];
          forcing[cell][x2] += f[1];
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& s,
        const SparseMatrixPolicy&) const
    {
      const auto x1 = s.at("X1"), x2 = s.at("X2");
      // ProcessSet invokes external Jacobian functions as (parameters, state, jacobian).
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          const double v1 = state[block][x1] - Q2_SHIFT[0], v2 = state[block][x2] - Q2_SHIFT[1];
          jacobian[block][x1][x1] -= -2.0 * v1 + 4.0 * v2 - 1.0;
          jacobian[block][x1][x2] -= 4.0 * v1 + 2.0;
          jacobian[block][x2][x1] -= 2.0 * v1 + 1.0;
          jacobian[block][x2][x2] -= 4.0 * v2;
        }
      };
    }
  };

  // ---- QDE-3 ----------------------------------------------------------------
  constexpr double Q3_SHIFT[3] = { 6.0, 1.0, 1.0 };
  constexpr double Q3_X0[3] = { -0.5, 0.3, 0.2 };  // unshifted IC
  constexpr double Q3_TEND = 2.0;

  inline void Qde3Rhs(const double x[3], double f[3])
  {
    f[0] = x[0] * x[0] + 7.0 * x[1] * x[1] - 4.0 * x[0] * x[1] + 5.0 * x[0];
    f[1] = -2.0 * x[1] * x[1] + x[0] * x[1] + 2.0 * x[0] * x[2] - 7.0 * x[1] * x[2] + 2.0 * x[1];
    f[2] = -x[1] * x[1] - 7.0 * x[2] * x[2] - 4.0 * x[1] * x[2] - x[2];
  }

  /// Exact shifted solution of QDE-3 via Bacsi-Kocsis eq. (17); M is diagonal.
  std::array<double, 3> Qde3Exact(double t)
  {
    const double E[3] = { std::exp(t), std::exp(-2.0 * t), std::exp(-5.0 * t) };
    const double Minv_w[3] = { 7.0, -1.0, 0.2 };
    double yp[3], ex[3];
    for (int i = 0; i < 3; ++i)
    {
      yp[i] = (E[i] - 1.0) * Minv_w[i];
      ex[i] = E[i] * Q3_X0[i];
    }
    // x^T B x = x1*x3 + x2^2 for B = [[0,0,1/2],[0,1,0],[1/2,0,0]]
    auto quad = [](const double u[3], const double v[3])
    { return 0.5 * (u[0] * v[2] + u[2] * v[0]) + u[1] * v[1]; };
    const double lambda = 4.0;
    const double b0 = quad(Q3_X0, Q3_X0);
    const double den = std::exp(-lambda * t) + 2.0 * quad(yp, ex) + b0 * quad(yp, yp);
    return { (ex[0] + b0 * yp[0]) / den + Q3_SHIFT[0],
             (ex[1] + b0 * yp[1]) / den + Q3_SHIFT[1],
             (ex[2] + b0 * yp[2]) / den + Q3_SHIFT[2] };
  }

  /// QDE-3 in shifted coordinates as a MICM external model.
  class Qde3Model
  {
   public:
    std::set<std::string> SpeciesUsed() const
    {
      return { "X1", "X2", "X3" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto x1 = s.at("X1"), x2 = s.at("X2"), x3 = s.at("X3");
      return { { x1, x1 }, { x1, x2 }, { x2, x1 }, { x2, x2 }, { x2, x3 }, { x3, x2 }, { x3, x3 } };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> UpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) { };
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto x1 = s.at("X1"), x2 = s.at("X2"), x3 = s.at("X3");
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          const double x[3] = { state[cell][x1] - Q3_SHIFT[0],
                                state[cell][x2] - Q3_SHIFT[1],
                                state[cell][x3] - Q3_SHIFT[2] };
          double f[3];
          Qde3Rhs(x, f);
          forcing[cell][x1] += f[0];
          forcing[cell][x2] += f[1];
          forcing[cell][x3] += f[2];
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& s,
        const SparseMatrixPolicy&) const
    {
      const auto x1 = s.at("X1"), x2 = s.at("X2"), x3 = s.at("X3");
      // ProcessSet invokes external Jacobian functions as (parameters, state, jacobian).
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, SparseMatrixPolicy& jacobian)
      {
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          const double v1 = state[block][x1] - Q3_SHIFT[0];
          const double v2 = state[block][x2] - Q3_SHIFT[1];
          const double v3 = state[block][x3] - Q3_SHIFT[2];
          jacobian[block][x1][x1] -= 2.0 * v1 - 4.0 * v2 + 5.0;
          jacobian[block][x1][x2] -= 14.0 * v2 - 4.0 * v1;
          jacobian[block][x2][x1] -= v2 + 2.0 * v3;
          jacobian[block][x2][x2] -= v1 - 4.0 * v2 - 7.0 * v3 + 2.0;
          jacobian[block][x2][x3] -= 2.0 * v1 - 7.0 * v2;
          jacobian[block][x3][x2] -= -2.0 * v2 - 4.0 * v3;
          jacobian[block][x3][x3] -= -4.0 * v2 - 14.0 * v3 - 1.0;
        }
      };
    }
  };

  // ---- harness --------------------------------------------------------------
  struct RunOutcome
  {
    std::uint64_t accepted = 0;
    std::uint64_t rejected = 0;
    bool converged = true;
    double max_rel_err = 0.0;
    double final_rel_err = 0.0;
  };

  /// Verify the closed form differentially: d/dt of the shifted exact solution
  /// must reproduce the shifted RHS. Catches transcription errors in either.
  bool SelfCheck()
  {
    const double eps = 1.0e-6;
    double worst = 0.0;
    for (double t : { 0.3, 1.1, 2.7 })
    {
      const auto xp = Qde2Exact(t + eps), xm = Qde2Exact(t - eps), xc = Qde2Exact(t);
      const double x[2] = { xc[0] - Q2_SHIFT[0], xc[1] - Q2_SHIFT[1] };
      double f[2];
      Qde2Rhs(x, f);
      for (int i = 0; i < 2; ++i)
        worst = std::max(worst, std::abs((xp[i] - xm[i]) / (2.0 * eps) - f[i]) / (std::abs(f[i]) + 1.0));
    }
    for (double t : { 0.3, 1.1, 1.8 })
    {
      const auto xp = Qde3Exact(t + eps), xm = Qde3Exact(t - eps), xc = Qde3Exact(t);
      const double x[3] = { xc[0] - Q3_SHIFT[0], xc[1] - Q3_SHIFT[1], xc[2] - Q3_SHIFT[2] };
      double f[3];
      Qde3Rhs(x, f);
      for (int i = 0; i < 3; ++i)
        worst = std::max(worst, std::abs((xp[i] - xm[i]) / (2.0 * eps) - f[i]) / (std::abs(f[i]) + 1.0));
    }
    std::cout << "exact-solution self-check residual: " << worst << "\n";
    return worst < 1.0e-5;
  }

  template<typename Solver, typename ExactFn>
  RunOutcome RunAgainstExact(
      Solver& solver,
      const std::vector<std::pair<std::string, double>>& initial_values,
      double rtol,
      double atol,
      double t_end,
      int n_output,
      ExactFn exact)
  {
    auto state = solver.GetState(1);
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, atol));
    for (const auto& [name, value] : initial_values)
      state.variables_[0][state.variable_map_.at(name)] = value;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;
    solver.UpdateStateParameters(state);

    RunOutcome out;
    double current = 0.0;
    for (int i_out = 1; i_out <= n_output; ++i_out)
    {
      const double t_out = t_end * i_out / n_output;
      double done = 0.0;
      const double dt = t_out - current;
      while (done < dt)
      {
        auto result = solver.Solve(dt - done, state);
        if (result.state_ != micm::SolverState::Converged && result.state_ != micm::SolverState::ConvergenceExceededMaxSteps)
        {
          out.converged = false;
          return out;
        }
        if (result.stats_.final_time_ <= 0.0)
        {
          out.converged = false;
          return out;
        }
        out.accepted += result.stats_.accepted_;
        out.rejected += result.stats_.rejected_;
        done += result.stats_.final_time_;
      }
      current = t_out;
      const auto xe = exact(t_out);
      double err = 0.0;
      for (std::size_t k = 0; k < xe.size(); ++k)
      {
        const double got = state.variables_[0][state.variable_map_.at("X" + std::to_string(k + 1))];
        err = std::max(err, std::abs(got - xe[k]) / (std::abs(xe[k]) + 1e-30));
      }
      out.max_rel_err = std::max(out.max_rel_err, err);
      if (i_out == n_output)
        out.final_rel_err = err;
    }
    return out;
  }

  template<typename ModelType, typename ExactFn>
  void RunProblem(
      const std::string& problem,
      const ModelType& model,
      const std::vector<std::string>& species,
      const std::vector<std::pair<std::string, double>>& init,
      double t_end,
      ExactFn exact,
      std::ofstream& tol_csv,
      std::ofstream& fix_csv)
  {
    const std::vector<double> rtols = { 1e-4, 1e-6, 1e-8, 1e-10 };
    for (double rtol : rtols)
    {
      std::vector<micm::PhaseSpecies> phase_species;
      for (const auto& name : species)
        phase_species.push_back(micm::Species(name));
      const micm::Phase gas{ "gas", phase_species };
      auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
      auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                        .SetSystem(micm::System(gas))
                        .AddExternalModel(model)
                        .SetReorderState(false)
                        .Build();
      const auto r = RunAgainstExact(solver, init, rtol, 1e-14, t_end, 20, exact);
      tol_csv << problem << ",four_stage," << rtol << ',' << r.accepted << ',' << r.rejected << ','
              << (r.converged ? 1 : 0) << ',' << r.max_rel_err << '\n';
      std::cout << problem << " four_stage rtol=" << rtol << "  accepted=" << r.accepted
                << "  max_rel_err=" << r.max_rel_err << (r.converged ? "" : "  [FAILED]") << "\n";
    }

    const std::vector<std::pair<std::string, std::function<micm::RosenbrockSolverParameters()>>> methods = {
      { "four_stage", []() { return micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters(); } },
      { "six_stage", []() { return micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters(); } },
      { "rodas4p", []() { return micm::RosenbrockSolverParameters::Rodas4PDifferentialAlgebraicRosenbrockParameters(); } },
    };
    for (const auto& [method_name, make_options] : methods)
    {
      for (int k = 0; k <= 7; ++k)
      {
        const double H = (t_end / 10.0) / std::pow(2.0, k);
        std::vector<micm::PhaseSpecies> phase_species;
        for (const auto& name : species)
          phase_species.push_back(micm::Species(name));
        const micm::Phase gas{ "gas", phase_species };
        auto options = make_options();
        options.h_start_ = H;
        options.h_min_ = H;
        options.h_max_ = H;
        auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                          .SetSystem(micm::System(gas))
                          .AddExternalModel(model)
                          .SetReorderState(false)
                          .Build();
        // Sample the fixed-step error at three independent horizons (each an
        // exact multiple of H), taking the max, so the order estimate is not
        // hostage to a sign change of the leading error term at one sample
        // time. Rows with rejected > 0 are not truly fixed-step (the
        // controller sub-stepped) and are excluded from order fits.
        RunOutcome r_full;
        double err_max = 0.0;
        bool converged = true;
        for (double frac : { 0.6, 0.8, 1.0 })
        {
          const auto r = RunAgainstExact(solver, init, 1e-6, 1e-14, frac * t_end, 1, exact);
          err_max = std::max(err_max, r.final_rel_err);
          converged = converged && r.converged;
          if (frac == 1.0)
            r_full = r;
        }
        fix_csv << problem << ',' << method_name << ',' << H << ',' << r_full.accepted << ',' << r_full.rejected << ','
                << (converged ? 1 : 0) << ',' << err_max << ',' << r_full.final_rel_err << '\n';
        std::cout << problem << " " << method_name << " H=" << H << "  err_max=" << err_max
                  << "  rejected=" << r_full.rejected << (converged ? "" : "  [FAILED]") << "\n";
      }
    }
  }
}  // namespace

int main()
{
  if (!SelfCheck())
  {
    std::cout << "closed-form self-check FAILED\n";
    return 1;
  }

  std::ofstream tol_csv("qde_exact.csv");
  tol_csv << "problem,method,rtol,accepted,rejected,converged,max_rel_err\n";
  std::ofstream fix_csv("qde_fixed_step.csv");
  fix_csv << "problem,method,H,steps,rejected,converged,err_max,err_final\n";

  RunProblem(
      "qde2",
      Qde2Model(),
      { "X1", "X2" },
      { { "X1", Q2_X0[0] + Q2_SHIFT[0] }, { "X2", Q2_X0[1] + Q2_SHIFT[1] } },
      Q2_TEND,
      [](double t) { return Qde2Exact(t); },
      tol_csv,
      fix_csv);

  RunProblem(
      "qde3",
      Qde3Model(),
      { "X1", "X2", "X3" },
      { { "X1", Q3_X0[0] + Q3_SHIFT[0] }, { "X2", Q3_X0[1] + Q3_SHIFT[1] }, { "X3", Q3_X0[2] + Q3_SHIFT[2] } },
      Q3_TEND,
      [](double t) { return Qde3Exact(t); },
      tol_csv,
      fix_csv);

  std::cout << "wrote qde_exact.csv and qde_fixed_step.csv\n";
  return 0;
}
