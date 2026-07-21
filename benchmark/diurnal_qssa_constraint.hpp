// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Time-dependent variant of the tropospheric QSSA radical constraint: the
// three photolysis rates are read per evaluation from the chemistry's custom
// rate parameters ("p1", "p2", "p3") instead of being baked in at
// construction, so the same solver instance follows a diurnal J(t) cycle set
// segment-by-segment by the driver. Thermal rate constants and the per-row
// residual scaling (computed from a reference — e.g. noon — radical state)
// remain fixed; scaling a complete constraint row by a constant leaves the
// solution unchanged.
#pragma once

#include "tropospheric_qssa_constraint.hpp"
#include "tropospheric_system.hpp"

#include <micm/system/conditions.hpp>

#include <functional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace tropospheric
{
  class DiurnalQssaRadicalConstraint
  {
   public:
    /// `reference_rates` supplies the thermal ks and, via `reference`, the
    /// per-row scaling; its j fields are ignored at evaluation time.
    DiurnalQssaRadicalConstraint(RadicalRates reference_rates, const RefConc& reference)
        : r_(reference_rates)
    {
      const double L_O1D = r_.k_O1D_N2 * reference.N2 + r_.k_O1D_O2 * reference.O2 + r_.k_O1D_H2O * reference.H2O;
      const double L_O = r_.k_O_O2_M * reference.O2 * reference.M + r_.k_O_O3 * reference.O3;
      const double L_OH = r_.k_OH_CO * reference.CO + r_.k_OH_CH4 * reference.CH4 + r_.k_OH_NO2 * reference.NO2;
      const double L_HO2 = r_.k_HO2_NO * reference.NO + 2.0 * r_.k_HO2_HO2 * reference.HO2;
      s_O1D_ = 1.0 / (L_O1D * reference.O1D);
      s_O_ = 1.0 / (L_O * reference.O);
      s_OH_ = 1.0 / (L_OH * reference.OH);
      s_HO2_ = 1.0 / (L_HO2 * reference.HO2);
    }

    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      return { "O1D", "O", "OH", "HO2" };
    }

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
        { O1D, O1D }, { O1D, O3 }, { O1D, N2 },  { O1D, O2 }, { O1D, H2O },
        { O, O },     { O, NO2 },  { O, O3 },    { O, O1D },  { O, N2 },    { O, O2 },  { O, M },
        { OH, OH },   { OH, O1D }, { OH, H2O },  { OH, HO2 }, { OH, NO },   { OH, CO }, { OH, CH4 }, { OH, NO2 },
        { HO2, HO2 }, { HO2, OH }, { HO2, CO },  { HO2, CH4 }, { HO2, NO }
      };
    }

    std::set<std::string> ConstraintStateParameterNames() const
    {
      return {};
    }

    template<typename DenseMatrixPolicy>
    std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)> ConstraintUpdateStateParametersFunction(
        const std::unordered_map<std::string, std::size_t>&) const
    {
      return [](const std::vector<micm::Conditions>&, DenseMatrixPolicy&) {};
    }

    // Residuals with per-cell photolysis read from the parameters matrix.
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::unordered_map<std::string, std::size_t>& parameter_indices,
        const std::unordered_map<std::string, std::size_t>& v) const
    {
      const auto iO1D = v.at("O1D"), iO = v.at("O"), iOH = v.at("OH"), iHO2 = v.at("HO2");
      const auto iO3 = v.at("O3"), iNO = v.at("NO"), iNO2 = v.at("NO2"), iCO = v.at("CO"), iCH4 = v.at("CH4");
      const auto iN2 = v.at("N2"), iO2 = v.at("O2"), iM = v.at("M"), iH2O = v.at("H2O");
      const auto ip1 = parameter_indices.at("p1"), ip2 = parameter_indices.at("p2"), ip3 = parameter_indices.at("p3");
      RadicalRates r = r_;
      const double s_O1D = s_O1D_, s_O = s_O_, s_OH = s_OH_, s_HO2 = s_HO2_;
      return [=](const DenseMatrixPolicy& y, const DenseMatrixPolicy& p, DenseMatrixPolicy& f)
      {
        for (std::size_t i = 0; i < y.NumRows(); ++i)
        {
          const double j1 = p[i][ip1], j2 = p[i][ip2], j3 = p[i][ip3];
          const double O1D = y[i][iO1D], O = y[i][iO], OH = y[i][iOH], HO2 = y[i][iHO2];
          const double O3 = y[i][iO3], NO = y[i][iNO], NO2 = y[i][iNO2], CO = y[i][iCO], CH4 = y[i][iCH4];
          const double N2 = y[i][iN2], O2 = y[i][iO2], M = y[i][iM], H2O = y[i][iH2O];

          f[i][iO1D] = s_O1D * (j2 * O3 - (r.k_O1D_N2 * N2 + r.k_O1D_O2 * O2 + r.k_O1D_H2O * H2O) * O1D);
          f[i][iO] = s_O * (j1 * NO2 + j3 * O3 + (r.k_O1D_N2 * N2 + r.k_O1D_O2 * O2) * O1D -
                            (r.k_O_O2_M * O2 * M + r.k_O_O3 * O3) * O);
          f[i][iOH] = s_OH * (2.0 * r.k_O1D_H2O * H2O * O1D + r.k_HO2_NO * NO * HO2 -
                              (r.k_OH_CO * CO + r.k_OH_CH4 * CH4 + r.k_OH_NO2 * NO2) * OH);
          f[i][iHO2] =
              s_HO2 * ((r.k_OH_CO * CO + r.k_OH_CH4 * CH4) * OH - r.k_HO2_NO * NO * HO2 - 2.0 * r.k_HO2_HO2 * HO2 * HO2);
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>& parameter_indices,
        const std::unordered_map<std::string, std::size_t>& v,
        const SparseMatrixPolicy&) const
    {
      const auto iO1D = v.at("O1D"), iO = v.at("O"), iOH = v.at("OH"), iHO2 = v.at("HO2");
      const auto iO3 = v.at("O3"), iNO = v.at("NO"), iNO2 = v.at("NO2"), iCO = v.at("CO"), iCH4 = v.at("CH4");
      const auto iN2 = v.at("N2"), iO2 = v.at("O2"), iM = v.at("M"), iH2O = v.at("H2O");
      const auto ip1 = parameter_indices.at("p1"), ip2 = parameter_indices.at("p2"), ip3 = parameter_indices.at("p3");
      RadicalRates r = r_;
      const double s_O1D = s_O1D_, s_O = s_O_, s_OH = s_OH_, s_HO2 = s_HO2_;
      return [=](const DenseMatrixPolicy& y, const DenseMatrixPolicy& p, SparseMatrixPolicy& jac)
      {
        for (std::size_t i = 0; i < jac.NumberOfBlocks(); ++i)
        {
          const double j1 = p[i][ip1], j2 = p[i][ip2], j3 = p[i][ip3];
          const double O1D = y[i][iO1D], O = y[i][iO], OH = y[i][iOH], HO2 = y[i][iHO2];
          const double O3 = y[i][iO3], NO = y[i][iNO], CO = y[i][iCO], CH4 = y[i][iCH4];
          const double N2 = y[i][iN2], O2 = y[i][iO2], M = y[i][iM], H2O = y[i][iH2O];

          jac[i][iO1D][iO1D] -= s_O1D * -(r.k_O1D_N2 * N2 + r.k_O1D_O2 * O2 + r.k_O1D_H2O * H2O);
          jac[i][iO1D][iO3] -= s_O1D * j2;
          jac[i][iO1D][iN2] -= s_O1D * -r.k_O1D_N2 * O1D;
          jac[i][iO1D][iO2] -= s_O1D * -r.k_O1D_O2 * O1D;
          jac[i][iO1D][iH2O] -= s_O1D * -r.k_O1D_H2O * O1D;

          jac[i][iO][iO] -= s_O * -(r.k_O_O2_M * O2 * M + r.k_O_O3 * O3);
          jac[i][iO][iNO2] -= s_O * j1;
          jac[i][iO][iO3] -= s_O * (j3 - r.k_O_O3 * O);
          jac[i][iO][iO1D] -= s_O * (r.k_O1D_N2 * N2 + r.k_O1D_O2 * O2);
          jac[i][iO][iN2] -= s_O * r.k_O1D_N2 * O1D;
          jac[i][iO][iO2] -= s_O * (r.k_O1D_O2 * O1D - r.k_O_O2_M * O * M);
          jac[i][iO][iM] -= s_O * -r.k_O_O2_M * O * O2;

          jac[i][iOH][iOH] -= s_OH * -(r.k_OH_CO * CO + r.k_OH_CH4 * CH4 + r.k_OH_NO2 * y[i][iNO2]);
          jac[i][iOH][iO1D] -= s_OH * 2.0 * r.k_O1D_H2O * H2O;
          jac[i][iOH][iH2O] -= s_OH * 2.0 * r.k_O1D_H2O * O1D;
          jac[i][iOH][iHO2] -= s_OH * r.k_HO2_NO * NO;
          jac[i][iOH][iNO] -= s_OH * r.k_HO2_NO * HO2;
          jac[i][iOH][iCO] -= s_OH * -r.k_OH_CO * OH;
          jac[i][iOH][iCH4] -= s_OH * -r.k_OH_CH4 * OH;
          jac[i][iOH][iNO2] -= s_OH * -r.k_OH_NO2 * OH;

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
    double s_O1D_ = 1.0, s_O_ = 1.0, s_OH_ = 1.0, s_HO2_ = 1.0;
  };
}  // namespace tropospheric
