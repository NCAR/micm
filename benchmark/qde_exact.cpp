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

  // ---- QDE-S: stiff rotated family ------------------------------------------
  // Inverse-designed stiff member of the same solvable class:
  //   V = R diag(-kappa, -1, -0.1) R^T with R = R12(4/5,3/5) R13(4/5,3/5)
  //   B = sym(r2 r3^T) (the two SLOW eigenvectors, so x^T B x stays O(1))
  //   lambda = s2 + s3 = -1.1, w = (0.4, 0.8, -0.6), A_i from their eq. (13):
  //   A_i = w_i B - (Bw) e_i^T - e_i (Bw)^T.
  // The rotation makes every component fully two-way coupled; the slow
  // dynamics are kappa-independent by construction (kappa lives only in the
  // fast eigendirection), so kappa is a pure stiffness knob on one and the
  // same trajectory. The initial condition sits on the slow manifold
  // (y0 = r2 + 2 r3 in inverse coordinates => x0 = (r2 + 2 r3)/2, exact
  // decimals from the 3-4-5 rotations), so there is no initial layer.
  using Vec3 = std::array<double, 3>;
  using Mat3 = std::array<Vec3, 3>;

  constexpr double QS_R[3][3] = { { 0.64, -0.6, -0.48 }, { 0.48, 0.8, -0.36 }, { 0.6, 0.0, 0.8 } };
  constexpr double QS_LAMBDA = -1.1;
  constexpr double QS_W[3] = { 0.4, 0.8, -0.6 };
  constexpr double QS_SHIFT[3] = { 2.0, 1.0, 1.0 };
  constexpr double QS_X0[3] = { -0.78, 0.04, 0.8 };  // = (r2 + 2 r3)/2, on the slow manifold
  constexpr double QS_TEND = 5.0;

  const Mat3& QsB()
  {
    static const Mat3 B = []()
    {
      Mat3 b{};
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          b[i][j] = 0.5 * (QS_R[i][1] * QS_R[j][2] + QS_R[i][2] * QS_R[j][1]);
      return b;
    }();
    return B;
  }

  const std::array<Mat3, 3>& QsA()
  {
    static const std::array<Mat3, 3> A = []()
    {
      const Mat3& b = QsB();
      Vec3 bw{};
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          bw[i] += b[i][j] * QS_W[j];
      std::array<Mat3, 3> a{};
      for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c)
          {
            a[i][r][c] = QS_W[i] * b[r][c];
            if (c == i)
              a[i][r][c] -= bw[r];
            if (r == i)
              a[i][r][c] -= bw[c];
          }
      return a;
    }();
    return A;
  }

  Mat3 QsV(double kappa)
  {
    const double s[3] = { -kappa, -1.0, -0.1 };
    Mat3 V{};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          V[i][j] += QS_R[i][k] * s[k] * QS_R[j][k];
    return V;
  }

  /// Original-coordinate RHS of QDE-S at stiffness kappa.
  void QdeStiffRhs(double kappa, const double x[3], double f[3])
  {
    const Mat3 V = QsV(kappa);
    const auto& A = QsA();
    for (int i = 0; i < 3; ++i)
    {
      f[i] = 0.0;
      for (int r = 0; r < 3; ++r)
      {
        f[i] += V[i][r] * x[r];
        for (int c = 0; c < 3; ++c)
          f[i] += A[i][r][c] * x[r] * x[c];
      }
    }
  }

  /// Exact shifted solution of QDE-S via Bacsi-Kocsis eq. (17);
  /// e^{Mt} = R diag(e^{mu_k t}) R^T with mu_k = s_k - lambda.
  std::array<double, 3> QdeStiffExact(double kappa, double t)
  {
    const Mat3& B = QsB();
    const double mu[3] = { -kappa - QS_LAMBDA, -1.0 - QS_LAMBDA, -0.1 - QS_LAMBDA };
    double rtw[3] = { 0.0, 0.0, 0.0 }, minvw[3] = { 0.0, 0.0, 0.0 };
    for (int k = 0; k < 3; ++k)
      for (int i = 0; i < 3; ++i)
        rtw[k] += QS_R[i][k] * QS_W[i];
    for (int i = 0; i < 3; ++i)
      for (int k = 0; k < 3; ++k)
        minvw[i] += QS_R[i][k] * rtw[k] / mu[k];
    double E[3][3];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        E[i][j] = 0.0;
        for (int k = 0; k < 3; ++k)
          E[i][j] += QS_R[i][k] * std::exp(mu[k] * t) * QS_R[j][k];
      }
    double yp[3], ex[3];
    for (int i = 0; i < 3; ++i)
    {
      yp[i] = -minvw[i];
      ex[i] = 0.0;
      for (int j = 0; j < 3; ++j)
      {
        yp[i] += E[i][j] * minvw[j];
        ex[i] += E[i][j] * QS_X0[j];
      }
    }
    auto quad = [&B](const double u[3], const double v[3])
    {
      double q = 0.0;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          q += u[i] * B[i][j] * v[j];
      return q;
    };
    const double b0 = quad(QS_X0, QS_X0);
    const double den = std::exp(-QS_LAMBDA * t) + 2.0 * quad(yp, ex) + b0 * quad(yp, yp);
    return { (ex[0] + b0 * yp[0]) / den + QS_SHIFT[0],
             (ex[1] + b0 * yp[1]) / den + QS_SHIFT[1],
             (ex[2] + b0 * yp[2]) / den + QS_SHIFT[2] };
  }

  /// QDE-S in shifted coordinates as a MICM external model, at fixed kappa.
  class QdeStiffModel
  {
   public:
    explicit QdeStiffModel(double kappa)
        : kappa_(kappa)
    {
    }

    std::set<std::string> SpeciesUsed() const
    {
      return { "X1", "X2", "X3" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto x1 = s.at("X1"), x2 = s.at("X2"), x3 = s.at("X3");
      std::set<std::pair<std::size_t, std::size_t>> elements;
      for (const auto row : { x1, x2, x3 })
        for (const auto col : { x1, x2, x3 })
          elements.insert({ row, col });
      return elements;
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
      const double kappa = kappa_;
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, DenseMatrixPolicy& forcing)
      {
        for (std::size_t cell = 0; cell < state.NumRows(); ++cell)
        {
          const double x[3] = { state[cell][x1] - QS_SHIFT[0],
                                state[cell][x2] - QS_SHIFT[1],
                                state[cell][x3] - QS_SHIFT[2] };
          double f[3];
          QdeStiffRhs(kappa, x, f);
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
      const Mat3 V = QsV(kappa_);
      // ProcessSet invokes external Jacobian functions as (parameters, state, jacobian).
      return [=](const DenseMatrixPolicy&, const DenseMatrixPolicy& state, SparseMatrixPolicy& jacobian)
      {
        const std::size_t idx[3] = { x1, x2, x3 };
        const auto& A = QsA();
        for (std::size_t block = 0; block < jacobian.NumberOfBlocks(); ++block)
        {
          const double x[3] = { state[block][x1] - QS_SHIFT[0],
                                state[block][x2] - QS_SHIFT[1],
                                state[block][x3] - QS_SHIFT[2] };
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
              double d = V[i][j];
              for (int k = 0; k < 3; ++k)
                d += 2.0 * A[i][j][k] * x[k];
              jacobian[block][idx[i]][idx[j]] -= d;
            }
        }
      };
    }

   private:
    double kappa_;
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
    for (double kappa : { 1.0e2, 1.0e4 })
      for (double t : { 0.3, 1.7, 4.5 })
      {
        const auto xp = QdeStiffExact(kappa, t + eps), xm = QdeStiffExact(kappa, t - eps),
                   xc = QdeStiffExact(kappa, t);
        const double x[3] = { xc[0] - QS_SHIFT[0], xc[1] - QS_SHIFT[1], xc[2] - QS_SHIFT[2] };
        double f[3];
        QdeStiffRhs(kappa, x, f);
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

  /// Integrate at working tolerance and dump computed vs exact trajectories.
  template<typename Solver, typename ExactFn>
  void DumpTimeseries(
      Solver& solver,
      const std::vector<std::pair<std::string, double>>& initial_values,
      double rtol,
      double t_end,
      int n_output,
      ExactFn exact,
      const std::string& problem,
      std::ofstream& ts_csv)
  {
    auto state = solver.GetState(1);
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1e-14));
    for (const auto& [name, value] : initial_values)
      state.variables_[0][state.variable_map_.at(name)] = value;
    state.conditions_[0].temperature_ = 298.0;
    state.conditions_[0].pressure_ = 101325.0;
    solver.UpdateStateParameters(state);

    const auto x0 = exact(0.0);
    for (std::size_t k = 0; k < x0.size(); ++k)
      ts_csv << problem << ",0," << (k + 1) << ',' << state.variables_[0][state.variable_map_.at("X" + std::to_string(k + 1))]
             << ',' << x0[k] << '\n';
    double current = 0.0;
    for (int i_out = 1; i_out <= n_output; ++i_out)
    {
      const double t_out = t_end * i_out / n_output;
      double done = 0.0;
      const double dt = t_out - current;
      while (done < dt)
      {
        auto result = solver.Solve(dt - done, state);
        if (result.stats_.final_time_ <= 0.0)
          return;
        done += result.stats_.final_time_;
      }
      current = t_out;
      const auto xe = exact(t_out);
      for (std::size_t k = 0; k < xe.size(); ++k)
        ts_csv << problem << ',' << t_out << ',' << (k + 1) << ','
               << state.variables_[0][state.variable_map_.at("X" + std::to_string(k + 1))] << ',' << xe[k] << '\n';
    }
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

  // ---- QDE-S: stiff rotated family ----
  {
    const std::vector<std::string> species = { "X1", "X2", "X3" };
    const std::vector<std::pair<std::string, double>> init = { { "X1", QS_X0[0] + QS_SHIFT[0] },
                                                               { "X2", QS_X0[1] + QS_SHIFT[1] },
                                                               { "X3", QS_X0[2] + QS_SHIFT[2] } };
    auto build_solver = [&species](const QdeStiffModel& model, micm::RosenbrockSolverParameters options)
    {
      std::vector<micm::PhaseSpecies> phase_species;
      for (const auto& name : species)
        phase_species.push_back(micm::Species(name));
      const micm::Phase gas{ "gas", phase_species };
      return micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
          .SetSystem(micm::System(gas))
          .AddExternalModel(model)
          .SetReorderState(false)
          .Build();
    };

    // Adaptive tolerance sweep at the stiffest kappa.
    {
      const double kappa = 1.0e6;
      for (double rtol : { 1e-4, 1e-6, 1e-8, 1e-10 })
      {
        auto solver =
            build_solver(QdeStiffModel(kappa), micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
        const auto r = RunAgainstExact(
            solver, init, rtol, 1e-14, QS_TEND, 20, [kappa](double t) { return QdeStiffExact(kappa, t); });
        tol_csv << "qde_stiff,four_stage," << rtol << ',' << r.accepted << ',' << r.rejected << ','
                << (r.converged ? 1 : 0) << ',' << r.max_rel_err << '\n';
        std::cout << "qde_stiff four_stage rtol=" << rtol << "  accepted=" << r.accepted
                  << "  max_rel_err=" << r.max_rel_err << (r.converged ? "" : "  [FAILED]") << "\n";
      }
    }

    // Stiffness sweep: same trajectory, kappa spanning four decades.
    std::ofstream kappa_csv("qde_stiff_kappa.csv");
    kappa_csv << "kappa,rtol,accepted,rejected,converged,max_rel_err\n";
    for (double kappa : { 1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6 })
    {
      auto solver =
          build_solver(QdeStiffModel(kappa), micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
      const auto r = RunAgainstExact(
          solver, init, 1e-6, 1e-14, QS_TEND, 20, [kappa](double t) { return QdeStiffExact(kappa, t); });
      kappa_csv << kappa << ",1e-6," << r.accepted << ',' << r.rejected << ',' << (r.converged ? 1 : 0) << ','
                << r.max_rel_err << '\n';
      std::cout << "qde_stiff kappa=" << kappa << "  accepted=" << r.accepted << "  max_rel_err=" << r.max_rel_err
                << (r.converged ? "" : "  [FAILED]") << "\n";
    }

    // Fixed-step order study in the deep-stiff regime.
    const std::vector<std::pair<std::string, std::function<micm::RosenbrockSolverParameters()>>> methods = {
      { "four_stage", []() { return micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters(); } },
      { "six_stage", []() { return micm::RosenbrockSolverParameters::SixStageDifferentialAlgebraicRosenbrockParameters(); } },
      { "rodas4p", []() { return micm::RosenbrockSolverParameters::Rodas4PDifferentialAlgebraicRosenbrockParameters(); } },
    };
    for (double kappa : { 1.0e4, 1.0e6 })
    {
      const std::string problem = kappa == 1.0e4 ? "qde_stiff_1e4" : "qde_stiff_1e6";
      for (const auto& [method_name, make_options] : methods)
      {
        for (int k = 0; k <= 7; ++k)
        {
          const double H = (QS_TEND / 10.0) / std::pow(2.0, k);
          auto options = make_options();
          options.h_start_ = H;
          options.h_min_ = H;
          options.h_max_ = H;
          auto solver = build_solver(QdeStiffModel(kappa), options);
          RunOutcome r_full;
          double err_max = 0.0;
          bool converged = true;
          for (double frac : { 0.6, 0.8, 1.0 })
          {
            const auto r = RunAgainstExact(
                solver, init, 1e-6, 1e-14, frac * QS_TEND, 1, [kappa](double t) { return QdeStiffExact(kappa, t); });
            err_max = std::max(err_max, r.final_rel_err);
            converged = converged && r.converged;
            if (frac == 1.0)
              r_full = r;
          }
          fix_csv << problem << ',' << method_name << ',' << H << ',' << r_full.accepted << ',' << r_full.rejected
                  << ',' << (converged ? 1 : 0) << ',' << err_max << ',' << r_full.final_rel_err << '\n';
          std::cout << problem << " " << method_name << " H=" << H << "  err_max=" << err_max
                    << "  rejected=" << r_full.rejected << (converged ? "" : "  [FAILED]") << "\n";
        }
      }
    }
  }

  // ---- computed-vs-exact trajectories for the paper figures ----
  {
    std::ofstream ts_csv("qde_timeseries.csv");
    ts_csv << "problem,t,component,computed,exact\n";
    auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    {
      std::vector<micm::PhaseSpecies> ps = { micm::Species("X1"), micm::Species("X2") };
      const micm::Phase gas{ "gas", ps };
      auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                        .SetSystem(micm::System(gas))
                        .AddExternalModel(Qde2Model())
                        .SetReorderState(false)
                        .Build();
      DumpTimeseries(
          solver,
          { { "X1", Q2_X0[0] + Q2_SHIFT[0] }, { "X2", Q2_X0[1] + Q2_SHIFT[1] } },
          1e-6,
          Q2_TEND,
          40,
          [](double t) { return Qde2Exact(t); },
          "qde2",
          ts_csv);
    }
    {
      std::vector<micm::PhaseSpecies> ps = { micm::Species("X1"), micm::Species("X2"), micm::Species("X3") };
      const micm::Phase gas{ "gas", ps };
      auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                        .SetSystem(micm::System(gas))
                        .AddExternalModel(Qde3Model())
                        .SetReorderState(false)
                        .Build();
      DumpTimeseries(
          solver,
          { { "X1", Q3_X0[0] + Q3_SHIFT[0] }, { "X2", Q3_X0[1] + Q3_SHIFT[1] }, { "X3", Q3_X0[2] + Q3_SHIFT[2] } },
          1e-6,
          Q3_TEND,
          40,
          [](double t) { return Qde3Exact(t); },
          "qde3",
          ts_csv);
    }
    {
      const double kappa = 1.0e6;
      std::vector<micm::PhaseSpecies> ps = { micm::Species("X1"), micm::Species("X2"), micm::Species("X3") };
      const micm::Phase gas{ "gas", ps };
      auto solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                        .SetSystem(micm::System(gas))
                        .AddExternalModel(QdeStiffModel(kappa))
                        .SetReorderState(false)
                        .Build();
      DumpTimeseries(
          solver,
          { { "X1", QS_X0[0] + QS_SHIFT[0] }, { "X2", QS_X0[1] + QS_SHIFT[1] }, { "X3", QS_X0[2] + QS_SHIFT[2] } },
          1e-6,
          QS_TEND,
          40,
          [kappa](double t) { return QdeStiffExact(kappa, t); },
          "qde_stiff",
          ts_csv);
    }
  }

  std::cout << "wrote qde_exact.csv, qde_fixed_step.csv, qde_stiff_kappa.csv, and qde_timeseries.csv\n";
  return 0;
}
