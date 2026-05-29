// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// QSSA constraint for the tropospheric O3-NOx-HOx mechanism: replaces the ODE
// rows of the four fast radicals (O1D, O, OH, HO2) with algebraic equations
// "net chemical production = 0". The four constraints are COUPLED (OH<->HO2
// cycling, O1D->OH source), so this is a single external model owning four
// algebraic rows rather than four independent constraints.
//
// MICM's ProcessSet skips algebraic rows when accumulating chemistry forcing
// (process_set.hpp), so the residual here must reconstruct each radical's full
// net production analytically; chemistry never writes to these rows.
//
// Rate constants are evaluated once at the (fixed) box temperature/pressure
// from the SAME tropospheric:: Arrhenius parameters the chemistry uses, so the
// DAE and full ODE are provably the same chemistry. Photolysis j's are passed
// in (already scaled by the benchmark's solar factor) and set identically on
// the chemistry side via SetCustomRateParameter.
#pragma once

#include "tropospheric_system.hpp"

#include <micm/process/rate_constant/rate_constant_functions.hpp>
#include <micm/system/conditions.hpp>

#include <array>
#include <cmath>
#include <functional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace tropospheric
{
  /// Rate constants needed by the radical residuals, all evaluated at the box
  /// (T, P). j1=jNO2, j2=jO1D, j3=jO3P are the (scaled) photolysis rates.
  struct RadicalRates
  {
    double k_O1D_N2, k_O1D_O2, k_O1D_H2O;  // O1D loss
    double k_O_O2_M, k_O_O3;               // O loss
    double k_O3_NO;                        // (not in radical rows, kept for completeness/IC)
    double k_OH_CO, k_OH_CH4;              // OH loss to HO2
    double k_HO2_NO;                       // HO2->OH (and OH source)
    double k_OH_NO2;                       // OH termination
    double k_HO2_HO2;                      // HO2 termination
    double j1, j2, j3;                     // jNO2, jO1D, jO3P

    static RadicalRates At(double T, double P, double jNO2, double jO1D, double jO3P)
    {
      using micm::CalculateArrhenius;
      return RadicalRates{ .k_O1D_N2 = CalculateArrhenius(AR_O1D_N2, T, P),
                           .k_O1D_O2 = CalculateArrhenius(AR_O1D_O2, T, P),
                           .k_O1D_H2O = CalculateArrhenius(AR_O1D_H2O, T, P),
                           .k_O_O2_M = CalculateArrhenius(AR_O_O2_M, T, P),
                           .k_O_O3 = CalculateArrhenius(AR_O_O3, T, P),
                           .k_O3_NO = CalculateArrhenius(AR_O3_NO, T, P),
                           .k_OH_CO = CalculateArrhenius(AR_OH_CO, T, P),
                           .k_OH_CH4 = CalculateArrhenius(AR_OH_CH4, T, P),
                           .k_HO2_NO = CalculateArrhenius(AR_HO2_NO, T, P),
                           .k_OH_NO2 = CalculateArrhenius(AR_OH_NO2, T, P),
                           .k_HO2_HO2 = CalculateArrhenius(AR_HO2_HO2, T, P),
                           .j1 = jNO2,
                           .j2 = jO1D,
                           .j3 = jO3P };
    }
  };

  // Reference concentrations used only to set per-row residual scaling.
  struct RefConc
  {
    double O1D, O, OH, HO2, O3, NO, NO2, CO, CH4, N2, O2, M, H2O;
  };

  class QssaRadicalConstraint
  {
   public:
    // The solver's constraint-initialization Newton converges on an ABSOLUTE
    // residual tolerance (||G||_inf < 1e-10). Our radical net-production terms
    // are ~1e8 molec/cm^3/s, so floating-point cancellation floors |G| near
    // 1e-7 and the absolute test can never pass. We therefore scale each
    // algebraic row by s_R = 1/(L_R * R_ref) (L_R = loss frequency, R_ref =
    // reference concentration), turning the residual into the RELATIVE error of
    // R. Scaling residual AND Jacobian rows by the same constant leaves the
    // Newton/Rosenbrock solution unchanged; only the convergence metric becomes
    // relative. (Robertson did not need this — its residuals were O(1).)
    QssaRadicalConstraint(RadicalRates r, const RefConc& ref)
        : r_(r)
    {
      const double L_O1D = r.k_O1D_N2 * ref.N2 + r.k_O1D_O2 * ref.O2 + r.k_O1D_H2O * ref.H2O;
      const double L_O = r.k_O_O2_M * ref.O2 * ref.M + r.k_O_O3 * ref.O3;
      const double L_OH = r.k_OH_CO * ref.CO + r.k_OH_CH4 * ref.CH4 + r.k_OH_NO2 * ref.NO2;
      const double L_HO2 = r.k_HO2_NO * ref.NO + 2.0 * r.k_HO2_HO2 * ref.HO2;
      s_O1D_ = 1.0 / (L_O1D * ref.O1D);
      s_O_ = 1.0 / (L_O * ref.O);
      s_OH_ = 1.0 / (L_OH * ref.OH);
      s_HO2_ = 1.0 / (L_HO2 * ref.HO2);
    }

    // The four fast-radical rows become algebraic.
    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      return { "O1D", "O", "OH", "HO2" };
    }

    // Every species appearing in any radical residual.
    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      return { "O1D", "O", "OH", "HO2", "O3", "NO", "NO2", "CO", "CH4", "N2", "O2", "M", "H2O" };
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      const auto O1D = s.at("O1D"), O = s.at("O"), OH = s.at("OH"), HO2 = s.at("HO2");
      const auto O3 = s.at("O3"), NO = s.at("NO"), NO2 = s.at("NO2"), CO = s.at("CO"), CH4 = s.at("CH4");
      const auto N2 = s.at("N2"), O2 = s.at("O2"), M = s.at("M"), H2O = s.at("H2O");
      return {
        // O1D row depends on O1D, O3, N2, O2, H2O
        { O1D, O1D }, { O1D, O3 }, { O1D, N2 }, { O1D, O2 }, { O1D, H2O },
        // O row depends on O, NO2, O3, O1D, N2, O2, M
        { O, O }, { O, NO2 }, { O, O3 }, { O, O1D }, { O, N2 }, { O, O2 }, { O, M },
        // OH row depends on OH, O1D, H2O, HO2, NO, CO, CH4, NO2
        { OH, OH }, { OH, O1D }, { OH, H2O }, { OH, HO2 }, { OH, NO }, { OH, CO }, { OH, CH4 }, { OH, NO2 },
        // HO2 row depends on HO2, OH, CO, CH4, NO
        { HO2, HO2 }, { HO2, OH }, { HO2, CO }, { HO2, CH4 }, { HO2, NO }
      };
    }

    // No constraint-owned state parameters: all rates are baked in at construction.
    std::set<std::string> ConstraintStateParameterNames() const
    {
      return {};
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>
    ConstraintUpdateStateParametersFunction(const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) {};
    }

    // Residuals G_R = net chemical production of R (forcing rows overwritten).
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& v) const
    {
      const auto iO1D = v.at("O1D"), iO = v.at("O"), iOH = v.at("OH"), iHO2 = v.at("HO2");
      const auto iO3 = v.at("O3"), iNO = v.at("NO"), iNO2 = v.at("NO2"), iCO = v.at("CO"), iCH4 = v.at("CH4");
      const auto iN2 = v.at("N2"), iO2 = v.at("O2"), iM = v.at("M"), iH2O = v.at("H2O");
      RadicalRates r = r_;
      const double s_O1D = s_O1D_, s_O = s_O_, s_OH = s_OH_, s_HO2 = s_HO2_;
      return [=](const DenseMatrixPolicy& y, const DenseMatrixPolicy&, DenseMatrixPolicy& f)
      {
        for (std::size_t i = 0; i < y.NumRows(); ++i)
        {
          const double O1D = y[i][iO1D], O = y[i][iO], OH = y[i][iOH], HO2 = y[i][iHO2];
          const double O3 = y[i][iO3], NO = y[i][iNO], NO2 = y[i][iNO2], CO = y[i][iCO], CH4 = y[i][iCH4];
          const double N2 = y[i][iN2], O2 = y[i][iO2], M = y[i][iM], H2O = y[i][iH2O];

          f[i][iO1D] = s_O1D * (r.j2 * O3 - (r.k_O1D_N2 * N2 + r.k_O1D_O2 * O2 + r.k_O1D_H2O * H2O) * O1D);
          f[i][iO] = s_O * (r.j1 * NO2 + r.j3 * O3 + (r.k_O1D_N2 * N2 + r.k_O1D_O2 * O2) * O1D -
                            (r.k_O_O2_M * O2 * M + r.k_O_O3 * O3) * O);
          f[i][iOH] = s_OH * (2.0 * r.k_O1D_H2O * H2O * O1D + r.k_HO2_NO * NO * HO2 -
                              (r.k_OH_CO * CO + r.k_OH_CH4 * CH4 + r.k_OH_NO2 * NO2) * OH);
          f[i][iHO2] =
              s_HO2 * ((r.k_OH_CO * CO + r.k_OH_CH4 * CH4) * OH - r.k_HO2_NO * NO * HO2 - 2.0 * r.k_HO2_HO2 * HO2 * HO2);
        }
      };
    }

    // Solver subtracts dG/dy: jac[row][col] -= dG_row/d(col).
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& v,
        const SparseMatrixPolicy&) const
    {
      const auto iO1D = v.at("O1D"), iO = v.at("O"), iOH = v.at("OH"), iHO2 = v.at("HO2");
      const auto iO3 = v.at("O3"), iNO = v.at("NO"), iNO2 = v.at("NO2"), iCO = v.at("CO"), iCH4 = v.at("CH4");
      const auto iN2 = v.at("N2"), iO2 = v.at("O2"), iM = v.at("M"), iH2O = v.at("H2O");
      RadicalRates r = r_;
      const double s_O1D = s_O1D_, s_O = s_O_, s_OH = s_OH_, s_HO2 = s_HO2_;
      return [=](const DenseMatrixPolicy& y, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
      {
        for (std::size_t i = 0; i < jac.NumberOfBlocks(); ++i)
        {
          const double O1D = y[i][iO1D], O = y[i][iO], OH = y[i][iOH], HO2 = y[i][iHO2];
          const double O3 = y[i][iO3], NO = y[i][iNO], NO2 = y[i][iNO2], CO = y[i][iCO], CH4 = y[i][iCH4];
          const double N2 = y[i][iN2], O2 = y[i][iO2], M = y[i][iM], H2O = y[i][iH2O];

          // G_O1D (row scaled by s_O1D)
          jac[i][iO1D][iO1D] -= s_O1D * -(r.k_O1D_N2 * N2 + r.k_O1D_O2 * O2 + r.k_O1D_H2O * H2O);
          jac[i][iO1D][iO3] -= s_O1D * r.j2;
          jac[i][iO1D][iN2] -= s_O1D * -r.k_O1D_N2 * O1D;
          jac[i][iO1D][iO2] -= s_O1D * -r.k_O1D_O2 * O1D;
          jac[i][iO1D][iH2O] -= s_O1D * -r.k_O1D_H2O * O1D;

          // G_O (row scaled by s_O)
          jac[i][iO][iO] -= s_O * -(r.k_O_O2_M * O2 * M + r.k_O_O3 * O3);
          jac[i][iO][iNO2] -= s_O * r.j1;
          jac[i][iO][iO3] -= s_O * (r.j3 - r.k_O_O3 * O);
          jac[i][iO][iO1D] -= s_O * (r.k_O1D_N2 * N2 + r.k_O1D_O2 * O2);
          jac[i][iO][iN2] -= s_O * r.k_O1D_N2 * O1D;
          jac[i][iO][iO2] -= s_O * (r.k_O1D_O2 * O1D - r.k_O_O2_M * O * M);
          jac[i][iO][iM] -= s_O * -r.k_O_O2_M * O * O2;

          // G_OH (row scaled by s_OH)
          jac[i][iOH][iOH] -= s_OH * -(r.k_OH_CO * CO + r.k_OH_CH4 * CH4 + r.k_OH_NO2 * NO2);
          jac[i][iOH][iO1D] -= s_OH * 2.0 * r.k_O1D_H2O * H2O;
          jac[i][iOH][iH2O] -= s_OH * 2.0 * r.k_O1D_H2O * O1D;
          jac[i][iOH][iHO2] -= s_OH * r.k_HO2_NO * NO;
          jac[i][iOH][iNO] -= s_OH * r.k_HO2_NO * HO2;
          jac[i][iOH][iCO] -= s_OH * -r.k_OH_CO * OH;
          jac[i][iOH][iCH4] -= s_OH * -r.k_OH_CH4 * OH;
          jac[i][iOH][iNO2] -= s_OH * -r.k_OH_NO2 * OH;

          // G_HO2 (row scaled by s_HO2)
          jac[i][iHO2][iHO2] -= s_HO2 * -(r.k_HO2_NO * NO + 4.0 * r.k_HO2_HO2 * HO2);
          jac[i][iHO2][iOH] -= s_HO2 * (r.k_OH_CO * CO + r.k_OH_CH4 * CH4);
          jac[i][iHO2][iCO] -= s_HO2 * r.k_OH_CO * OH;
          jac[i][iHO2][iCH4] -= s_HO2 * r.k_OH_CH4 * OH;
          jac[i][iHO2][iNO] -= s_HO2 * -r.k_HO2_NO * HO2;
        }
      };
    }

   private:
    RadicalRates r_;
    double s_O1D_ = 1.0, s_O_ = 1.0, s_OH_ = 1.0, s_HO2_ = 1.0;  // per-row relative-residual scales
  };

  /// Project the four radicals onto the QSSA manifold (G = 0) for the given
  /// slow-species concentrations, so the DAE starts consistent. O1D and O are
  /// exact (linear given upstream); OH/HO2 are coupled and solved by damped
  /// fixed-point iteration (HOx source = NO2 termination + HO2 self-reaction).
  struct ConsistentRadicals
  {
    double O1D, O, OH, HO2;
  };

  inline ConsistentRadicals
  ProjectRadicals(const RadicalRates& r, double O3, double NO, double NO2, double CO, double CH4, double N2, double O2, double M, double H2O)
  {
    ConsistentRadicals c{};
    c.O1D = (r.j2 * O3) / (r.k_O1D_N2 * N2 + r.k_O1D_O2 * O2 + r.k_O1D_H2O * H2O);
    c.O = (r.j1 * NO2 + r.j3 * O3 + (r.k_O1D_N2 * N2 + r.k_O1D_O2 * O2) * c.O1D) / (r.k_O_O2_M * O2 * M + r.k_O_O3 * O3);

    const double oh_loss = r.k_OH_CO * CO + r.k_OH_CH4 * CH4 + r.k_OH_NO2 * NO2;  // OH first-order loss
    const double oh_to_ho2 = r.k_OH_CO * CO + r.k_OH_CH4 * CH4;                   // OH->HO2 channel
    const double oh_src = 2.0 * r.k_O1D_H2O * H2O * c.O1D;                        // primary OH source
    double OH = oh_src / oh_loss;  // ignoring HO2->OH recycling (lower bound start)
    double HO2 = (oh_to_ho2 * OH) / (r.k_HO2_NO * NO + 1.0e-30);
    for (int it = 0; it < 500; ++it)
    {
      // HO2 from G_HO2=0: 2*k_HO2_HO2*HO2^2 + k_HO2_NO*NO*HO2 - oh_to_ho2*OH = 0 (positive root)
      const double a = 2.0 * r.k_HO2_HO2, b = r.k_HO2_NO * NO, cc = oh_to_ho2 * OH;
      const double HO2_new = (2.0 * cc) / (b + std::sqrt(b * b + 4.0 * a * cc));
      // OH from G_OH=0 with current HO2.
      const double OH_new = (oh_src + r.k_HO2_NO * NO * HO2_new) / oh_loss;
      OH = 0.5 * OH + 0.5 * OH_new;  // damping
      HO2 = 0.5 * HO2 + 0.5 * HO2_new;
    }
    c.OH = OH;
    c.HO2 = HO2;
    return c;
  }
}  // namespace tropospheric
