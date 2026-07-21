// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// MOZART-TS1 mechanism scale-up: the full imported TS1 gas-phase network
// (benchmark/data/ts1_mechanism.txt, built by benchmark/import_ts1.py from
// CAM's TS1-fullVBS chem_mech.in) integrated as a full ODE and as a QSSA-DAE
// with the classic fast-radical family {O1D, O, OH, HO2} replaced by a
// generic algebraic constraint computed from the reaction table itself.
//
// Conditions are a fixed polluted-boundary-layer box with steady clear-sky
// midday photolysis (major channels assigned representative values, the rest
// a small default), matching the steady-J design of the smaller tropospheric
// benchmark. Accuracy is measured on slow observables against a tight-
// tolerance full-ODE reference; wall-clock repetitions are interleaved.
#include <micm/CPU.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef MICM_BENCHMARK_DATA_DIR
#define MICM_BENCHMARK_DATA_DIR "benchmark/data"
#endif

namespace
{
  constexpr double T_BOX = 298.15;
  constexpr double P_BOX = 101325.0;
  constexpr double N_M = 2.5e19;

  // ---------------------------------------------------------------------------
  // Imported mechanism.
  // ---------------------------------------------------------------------------
  struct MechReaction
  {
    std::string kind;    // ARRH, TROE, PHOT
    std::string label;
    std::vector<double> params;
    std::vector<std::string> reactants;                    // with multiplicity
    std::vector<std::pair<double, std::string>> products;  // (coefficient, name)
  };

  struct Mechanism
  {
    std::vector<std::string> species;  // fixed reservoirs first
    std::set<std::string> fixed;
    std::vector<MechReaction> reactions;
  };

  Mechanism LoadMechanism(const std::string& path)
  {
    Mechanism mech;
    std::ifstream file(path);
    if (!file)
      throw std::runtime_error("cannot open " + path);
    std::string line;
    while (std::getline(file, line))
    {
      if (line.empty() || line[0] == '#')
        continue;
      std::istringstream ss(line);
      std::string tag;
      ss >> tag;
      if (tag == "S" || tag == "F")
      {
        std::string name;
        ss >> name;
        mech.species.push_back(name);
        if (tag == "F")
          mech.fixed.insert(name);
        continue;
      }
      // R KIND LABEL [params...] ; reactants ; products
      MechReaction reaction;
      ss >> reaction.kind >> reaction.label;
      std::string token;
      int section = 0;
      while (ss >> token)
      {
        if (token == ";")
        {
          ++section;
          continue;
        }
        if (section == 0)
        {
          reaction.params.push_back(std::stod(token));
        }
        else if (section == 1)
        {
          reaction.reactants.push_back(token);
        }
        else
        {
          const auto star = token.find('*');
          reaction.products.push_back({ std::stod(token.substr(0, star)), token.substr(star + 1) });
        }
      }
      mech.reactions.push_back(std::move(reaction));
    }
    return mech;
  }

  // Representative clear-sky midday photolysis rates (s^-1) for the major
  // channels; every other j gets a small default. A steady-J benchmark: the
  // values set the chemical regime, not the solver comparison.
  double PhotolysisRate(const std::string& label)
  {
    static const std::map<std::string, double> majors = {
      { "jo3_a", 3.0e-5 },   // O3 -> O1D
      { "jo3_b", 4.0e-4 },   // O3 -> O
      { "jno2", 8.0e-3 },    { "jno3_a", 1.5e-1 }, { "jno3_b", 2.0e-2 }, { "jhono", 1.5e-3 },
      { "jh2o2", 7.0e-6 },   { "jhno3", 5.0e-7 },  { "jho2no2_a", 3.5e-6 }, { "jho2no2_b", 1.5e-6 },
      { "jch2o_a", 3.0e-5 }, { "jch2o_b", 4.5e-5 }, { "jch3ooh", 5.0e-6 }, { "jn2o5_a", 3.0e-5 },
      { "jn2o5_b", 1.0e-6 }, { "jacet", 1.0e-6 },  { "jmvk", 2.0e-6 },   { "jmacr_a", 1.0e-6 },
      { "jglyoxal", 8.0e-5 }, { "jch3cho", 4.0e-6 },
    };
    // Unlisted channels (stratosphere-only photolysis like H2O/N2O, and
    // minor organic product channels) are OFF: a blanket default injects
    // enormous phantom radical sources at boundary-layer reservoir
    // concentrations (H2O + hv alone would dominate the O1D budget).
    const auto it = majors.find(label);
    return it != majors.end() ? it->second : 0.0;
  }

  // Initial number densities (cm^-3) for major species; everything else 0.
  double InitialValue(const std::string& name)
  {
    static const std::map<std::string, double> majors = {
      { "M", N_M },
      { "N2", 0.78 * N_M },
      { "O2", 0.21 * N_M },
      { "H2O", 0.01 * N_M },
      { "O3", 7.5e11 },   // 30 ppb
      { "NO", 2.5e10 },   // 1 ppb
      { "NO2", 2.5e10 },  // 1 ppb
      { "CO", 2.5e12 },   // 100 ppb
      { "CH4", 4.4e13 },  // 1.8 ppm
      { "H2", 1.25e13 },  // 500 ppb
      { "N2O", 8.2e12 },  // 330 ppb
      { "C2H6", 2.5e10 }, { "C3H8", 1.2e10 }, { "C2H4", 2.5e10 }, { "ISOP", 2.5e10 },
      { "CH2O", 5.0e10 }, { "CH3CHO", 2.5e10 }, { "CH3OH", 1.2e11 }, { "CH3COCH3", 2.5e10 },
      { "SO2", 2.5e10 },  { "DMS", 2.5e9 },   { "TOLUENE", 2.5e9 }, { "BENZENE", 2.5e9 },
    };
    const auto it = majors.find(name);
    return it != majors.end() ? it->second : 0.0;
  }

  /// Rate constant at the fixed box conditions.
  double RateAtBox(const MechReaction& reaction)
  {
    if (reaction.kind == "PHOT")
      return PhotolysisRate(reaction.label);
    if (reaction.kind == "ARRH")
    {
      micm::ArrheniusRateConstantParameters p;
      p.A_ = reaction.params[0];
      p.C_ = reaction.params.size() > 1 ? reaction.params[1] : 0.0;
      p.B_ = reaction.params.size() > 2 ? reaction.params[2] : 0.0;
      p.E_ = reaction.params.size() > 3 ? reaction.params[3] : 0.0;
      return micm::CalculateArrhenius(p, T_BOX, P_BOX);
    }
    micm::TroeRateConstantParameters p;
    p.k0_A_ = reaction.params[0];
    p.k0_B_ = reaction.params[1];
    p.kinf_A_ = reaction.params[2];
    p.kinf_B_ = reaction.params[3];
    p.Fc_ = reaction.params[4];
    return micm::CalculateTroe(p, T_BOX, N_M);
  }

  // ---------------------------------------------------------------------------
  // Generic QSSA constraint over the reaction table (fixed box conditions).
  // ---------------------------------------------------------------------------
  struct TableReaction
  {
    double k;
    std::vector<int> reactants;                    // species indices, with multiplicity
    std::vector<std::pair<double, int>> products;  // (coefficient, species index)
  };

  double NetStoichD(const TableReaction& r, int s)
  {
    double net = 0.0;
    for (const auto& [coef, p] : r.products)
      if (p == s)
        net += coef;
    for (int q : r.reactants)
      net -= (q == s);
    return net;
  }

  class Ts1QssaModel
  {
   public:
    Ts1QssaModel(
        std::vector<TableReaction> reactions,
        std::vector<std::string> species_names,
        std::vector<int> constrained,
        std::vector<double> row_scale)
        : reactions_(std::move(reactions)),
          species_names_(std::move(species_names)),
          constrained_(std::move(constrained)),
          row_scale_(std::move(row_scale))
    {
      // Reactions relevant to each constrained species.
      for (std::size_t c = 0; c < constrained_.size(); ++c)
      {
        std::vector<std::size_t> relevant;
        for (std::size_t j = 0; j < reactions_.size(); ++j)
        {
          if (NetStoichD(reactions_[j], constrained_[c]) != 0.0)
            relevant.push_back(j);
        }
        relevant_.push_back(std::move(relevant));
      }
    }

    std::set<std::string> ConstraintAlgebraicVariableNames() const
    {
      std::set<std::string> names;
      for (int c : constrained_)
        names.insert(species_names_[c]);
      return names;
    }

    std::set<std::string> ConstraintSpeciesDependencies() const
    {
      std::set<std::string> names;
      for (std::size_t c = 0; c < constrained_.size(); ++c)
      {
        names.insert(species_names_[constrained_[c]]);
        for (const std::size_t j : relevant_[c])
          for (int q : reactions_[j].reactants)
            names.insert(species_names_[q]);
      }
      return names;
    }

    std::set<std::pair<std::size_t, std::size_t>> NonZeroConstraintJacobianElements(
        const std::unordered_map<std::string, std::size_t>& s) const
    {
      std::set<std::pair<std::size_t, std::size_t>> elements;
      for (std::size_t c = 0; c < constrained_.size(); ++c)
      {
        const auto row = s.at(species_names_[constrained_[c]]);
        elements.insert({ row, row });
        for (const std::size_t j : relevant_[c])
          for (int q : reactions_[j].reactants)
            elements.insert({ row, s.at(species_names_[q]) });
      }
      return elements;
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

    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ConstraintResidualFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& v) const
    {
      auto self = *this;
      std::vector<std::size_t> index_of(species_names_.size());
      for (std::size_t i = 0; i < species_names_.size(); ++i)
        index_of[i] = v.at(species_names_[i]);
      return [self, index_of](const DenseMatrixPolicy& y, const DenseMatrixPolicy&, DenseMatrixPolicy& f)
      {
        for (std::size_t cell = 0; cell < y.NumRows(); ++cell)
        {
          for (std::size_t c = 0; c < self.constrained_.size(); ++c)
          {
            double net = 0.0;
            for (const std::size_t j : self.relevant_[c])
            {
              const auto& r = self.reactions_[j];
              double rate = r.k;
              for (int q : r.reactants)
                rate *= y[cell][index_of[q]];
              net += NetStoichD(r, self.constrained_[c]) * rate;
            }
            f[cell][index_of[self.constrained_[c]]] = self.row_scale_[c] * net;
          }
        }
      };
    }

    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> ConstraintJacobianFunction(
        const std::unordered_map<std::string, std::size_t>&,
        const std::unordered_map<std::string, std::size_t>& v,
        const SparseMatrixPolicy&) const
    {
      auto self = *this;
      std::vector<std::size_t> index_of(species_names_.size());
      for (std::size_t i = 0; i < species_names_.size(); ++i)
        index_of[i] = v.at(species_names_[i]);
      return [self, index_of](const DenseMatrixPolicy& y, const DenseMatrixPolicy&, SparseMatrixPolicy& jac)
      {
        for (std::size_t block = 0; block < jac.NumberOfBlocks(); ++block)
        {
          for (std::size_t c = 0; c < self.constrained_.size(); ++c)
          {
            const auto row = index_of[self.constrained_[c]];
            for (const std::size_t j : self.relevant_[c])
            {
              const auto& r = self.reactions_[j];
              const double stoich = NetStoichD(r, self.constrained_[c]);
              for (std::size_t which = 0; which < r.reactants.size(); ++which)
              {
                double derivative = r.k;
                for (std::size_t other = 0; other < r.reactants.size(); ++other)
                {
                  if (other != which)
                    derivative *= y[block][index_of[r.reactants[other]]];
                }
                jac[block][row][index_of[r.reactants[which]]] -= self.row_scale_[c] * stoich * derivative;
              }
            }
          }
        }
      };
    }

    /// Damped generic projection of the constrained species onto G = 0 with
    /// the other species held fixed (per-radical scalar Newton, cycled).
    void Project(std::vector<double>& values) const
    {
      for (int sweep = 0; sweep < 500; ++sweep)
      {
        for (std::size_t c = 0; c < constrained_.size(); ++c)
        {
          const int target = constrained_[c];
          double g = 0.0, dg = 0.0;
          for (const std::size_t j : relevant_[c])
          {
            const auto& r = reactions_[j];
            const double stoich = NetStoichD(r, target);
            double rate = r.k;
            for (int q : r.reactants)
              rate *= values[q];
            g += stoich * rate;
            for (std::size_t which = 0; which < r.reactants.size(); ++which)
            {
              if (r.reactants[which] != target)
                continue;
              double derivative = r.k;
              for (std::size_t other = 0; other < r.reactants.size(); ++other)
                if (other != which)
                  derivative *= values[r.reactants[other]];
              dg += stoich * derivative;
            }
          }
          if (dg != 0.0)
            values[target] = std::max(0.0, values[target] - 0.5 * g / dg);
        }
      }
    }

    std::vector<TableReaction> reactions_;
    std::vector<std::string> species_names_;
    std::vector<int> constrained_;
    std::vector<double> row_scale_;
    std::vector<std::vector<std::size_t>> relevant_;
  };

  // ---------------------------------------------------------------------------
  struct RunResult
  {
    std::uint64_t accepted = 0;
    std::uint64_t rejected = 0;
    bool converged = true;
    std::vector<double> o3, co, no2;
  };

  template<typename Solver>
  RunResult RunCase(
      Solver& solver,
      const Mechanism& mech,
      const std::vector<double>& initial,
      double rtol,
      const std::vector<double>& output_times)
  {
    auto state = solver.GetState(1);
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e2));
    for (const auto& reaction : mech.reactions)
    {
      if (reaction.kind == "PHOT")
        state.SetCustomRateParameter(reaction.label, PhotolysisRate(reaction.label));
    }
    auto& m = state.variable_map_;
    for (std::size_t i = 0; i < mech.species.size(); ++i)
      state.variables_[0][m.at(mech.species[i])] = initial[i];
    state.conditions_[0].temperature_ = T_BOX;
    state.conditions_[0].pressure_ = P_BOX;
    state.conditions_[0].air_density_ = N_M;
    solver.UpdateStateParameters(state);

    RunResult out;
    double current = 0.0;
    for (double t_out : output_times)
    {
      double done = 0.0;
      const double dt = t_out - current;
      while (done < dt)
      {
        auto result = solver.Solve(dt - done, state);
        if (result.state_ != micm::SolverState::Converged && result.state_ != micm::SolverState::ConvergenceExceededMaxSteps)
        {
          std::cout << "solver state: " << micm::SolverStateToString(result.state_) << " at t=" << t_out << "\n";
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
      out.o3.push_back(state.variables_[0][m.at("O3")]);
      out.co.push_back(state.variables_[0][m.at("CO")]);
      out.no2.push_back(state.variables_[0][m.at("NO2")]);
    }
    return out;
  }

  template<typename Solver>
  double TimedRun(
      Solver& solver,
      const Mechanism& mech,
      const std::vector<double>& initial,
      double rtol,
      const std::vector<double>& output_times)
  {
    auto state = solver.GetState(1);
    state.SetRelativeTolerance(rtol);
    state.SetAbsoluteTolerances(std::vector<double>(state.state_size_, 1.0e2));
    for (const auto& reaction : mech.reactions)
    {
      if (reaction.kind == "PHOT")
        state.SetCustomRateParameter(reaction.label, PhotolysisRate(reaction.label));
    }
    auto& m = state.variable_map_;
    for (std::size_t i = 0; i < mech.species.size(); ++i)
      state.variables_[0][m.at(mech.species[i])] = initial[i];
    state.conditions_[0].temperature_ = T_BOX;
    state.conditions_[0].pressure_ = P_BOX;
    state.conditions_[0].air_density_ = N_M;
    solver.UpdateStateParameters(state);
    const auto t0 = std::chrono::steady_clock::now();
    double current = 0.0;
    for (double t_out : output_times)
    {
      double done = 0.0;
      const double dt = t_out - current;
      while (done < dt)
      {
        auto result = solver.Solve(dt - done, state);
        if (result.state_ != micm::SolverState::Converged && result.state_ != micm::SolverState::ConvergenceExceededMaxSteps)
          break;
        if (result.stats_.final_time_ <= 0.0)
          break;
        done += result.stats_.final_time_;
      }
      current = t_out;
    }
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::micro>(t1 - t0).count();
  }

  double MaxRelError(const RunResult& a, const RunResult& b, const std::vector<double>& times, double t_skip)
  {
    double worst = 0.0;
    auto rel = [](double got, double r) { return std::abs(got - r) / (std::abs(r) + 1e-30); };
    for (std::size_t i = 0; i < times.size() && i < a.o3.size() && i < b.o3.size(); ++i)
    {
      if (times[i] < t_skip)
        continue;
      worst = std::max(worst, rel(a.o3[i], b.o3[i]));
      worst = std::max(worst, rel(a.co[i], b.co[i]));
      worst = std::max(worst, rel(a.no2[i], b.no2[i]));
    }
    return worst;
  }
}  // namespace

int main()
{
  const double rtol = 1e-6;
  const double t_skip = 10.0;
  const int reps = 11;
  std::vector<double> output_times;
  for (int i = 0; i < 13; ++i)
    output_times.push_back(std::pow(10.0, -2.0 + 7.0 * i / 12.0));  // 1e-2 .. 1e5 s

  auto mech = LoadMechanism(std::string(MICM_BENCHMARK_DATA_DIR) + "/ts1_mechanism.txt");
  std::cout << "TS1 (imported): " << mech.species.size() << " species, " << mech.reactions.size() << " reactions\n";

  // Species table and processes.
  std::unordered_map<std::string, int> species_index;
  std::vector<micm::PhaseSpecies> phase_species;
  for (std::size_t i = 0; i < mech.species.size(); ++i)
  {
    species_index[mech.species[i]] = static_cast<int>(i);
    phase_species.push_back(micm::Species(mech.species[i]));
  }
  micm::Phase gas{ "gas", phase_species };

  std::vector<micm::Process> processes;
  std::vector<TableReaction> table;
  for (const auto& reaction : mech.reactions)
  {
    std::vector<micm::Species> reactants;
    std::vector<micm::StoichSpecies> products;
    TableReaction table_reaction;
    table_reaction.k = RateAtBox(reaction);
    for (const auto& r : reaction.reactants)
    {
      reactants.push_back(micm::Species(r));
      table_reaction.reactants.push_back(species_index.at(r));
    }
    for (const auto& [coef, p] : reaction.products)
    {
      products.push_back({ micm::Species(p), coef });
      table_reaction.products.push_back({ coef, species_index.at(p) });
    }
    auto builder = micm::ChemicalReactionBuilder().SetReactants(reactants).SetProducts(products).SetPhase(gas);
    if (reaction.kind == "PHOT")
    {
      builder.SetRateConstant(micm::UserDefinedRateConstantParameters{ .label_ = reaction.label });
    }
    else if (reaction.kind == "ARRH")
    {
      micm::ArrheniusRateConstantParameters p;
      p.A_ = reaction.params[0];
      p.C_ = reaction.params.size() > 1 ? reaction.params[1] : 0.0;
      p.B_ = reaction.params.size() > 2 ? reaction.params[2] : 0.0;
      p.E_ = reaction.params.size() > 3 ? reaction.params[3] : 0.0;
      builder.SetRateConstant(p);
    }
    else
    {
      micm::TroeRateConstantParameters p;
      p.k0_A_ = reaction.params[0];
      p.k0_B_ = reaction.params[1];
      p.kinf_A_ = reaction.params[2];
      p.kinf_B_ = reaction.params[3];
      p.Fc_ = reaction.params[4];
      builder.SetRateConstant(p);
    }
    processes.push_back(builder.Build());
    table.push_back(std::move(table_reaction));
  }

  // Initial values, and the constrained radical family.
  std::vector<double> initial(mech.species.size());
  for (std::size_t i = 0; i < mech.species.size(); ++i)
    initial[i] = InitialValue(mech.species[i]);
  const std::vector<std::string> radical_names = { "O1D", "O", "OH", "HO2" };
  std::vector<int> constrained;
  for (const auto& name : radical_names)
    constrained.push_back(species_index.at(name));

  // Manifold projection and per-row scaling from the projected state.
  Ts1QssaModel prototype(table, mech.species, constrained, std::vector<double>(constrained.size(), 1.0));
  prototype.Project(initial);
  std::vector<double> row_scale;
  for (std::size_t c = 0; c < constrained.size(); ++c)
  {
    double loss = 0.0;
    for (const std::size_t j : prototype.relevant_[c])
    {
      const auto& r = prototype.reactions_[j];
      const double stoich = NetStoichD(r, constrained[c]);
      if (stoich >= 0.0)
        continue;
      double rate = r.k;
      for (int q : r.reactants)
        rate *= initial[q];
      loss += -stoich * rate;
    }
    row_scale.push_back(1.0 / std::max(loss, 1e-300));
  }
  Ts1QssaModel constraint(table, mech.species, constrained, row_scale);
  std::cout << "projected radicals: ";
  for (std::size_t c = 0; c < constrained.size(); ++c)
    std::cout << radical_names[c] << "=" << initial[constrained[c]] << " ";
  std::cout << "\n";

  auto options = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
  auto ode_solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                        .SetSystem(micm::System(gas))
                        .SetReactions(processes)
                        .SetReorderState(false)
                        .Build();
  auto dae_solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                        .SetSystem(micm::System(gas))
                        .SetReactions(processes)
                        .AddExternalModel(constraint)
                        .SetReorderState(false)
                        .Build();
  auto schur_options = options;
  schur_options.schur_reduction_ = true;
  auto schur_solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(schur_options)
                          .SetSystem(micm::System(gas))
                          .SetReactions(processes)
                          .AddExternalModel(constraint)
                          .SetReorderState(false)
                          .Build();

  auto reference = RunCase(ode_solver, mech, initial, 1e-9, output_times);
  auto ode = RunCase(ode_solver, mech, initial, rtol, output_times);
  auto dae = RunCase(dae_solver, mech, initial, rtol, output_times);
  auto dae_schur = RunCase(schur_solver, mech, initial, rtol, output_times);
  if (!reference.converged || !ode.converged || !dae.converged || !dae_schur.converged)
  {
    std::cout << "FAILED: reference=" << reference.converged << " ode=" << ode.converged << " dae=" << dae.converged
              << " schur=" << dae_schur.converged << "\n";
    return 1;
  }

  std::vector<double> ode_samples, dae_samples, schur_samples;
  for (int rep = 0; rep < reps; ++rep)
  {
    ode_samples.push_back(TimedRun(ode_solver, mech, initial, rtol, output_times));
    dae_samples.push_back(TimedRun(dae_solver, mech, initial, rtol, output_times));
    schur_samples.push_back(TimedRun(schur_solver, mech, initial, rtol, output_times));
  }
  auto median = [](std::vector<double>& v)
  {
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
  };
  const double ode_us = median(ode_samples);
  const double dae_us = median(dae_samples);
  const double schur_us = median(schur_samples);

  std::ofstream csv("ts1_dae.csv");
  csv << "method,rtol,accepted,rejected,wallclock_median_us,max_rel_err\n";
  csv.precision(12);
  const double ode_err = MaxRelError(ode, reference, output_times, t_skip);
  const double dae_err = MaxRelError(dae, reference, output_times, t_skip);
  const double schur_err = MaxRelError(dae_schur, reference, output_times, t_skip);
  csv << "full_ode," << rtol << ',' << ode.accepted << ',' << ode.rejected << ',' << ode_us << ',' << ode_err << '\n';
  csv << "dae_qssa," << rtol << ',' << dae.accepted << ',' << dae.rejected << ',' << dae_us << ',' << dae_err << '\n';
  csv << "dae_qssa_schur," << rtol << ',' << dae_schur.accepted << ',' << dae_schur.rejected << ',' << schur_us << ','
      << schur_err << '\n';

  std::cout << std::scientific << std::setprecision(3);
  std::cout << "ODE:       acc=" << ode.accepted << " rej=" << ode.rejected << " us=" << ode_us << " err=" << ode_err
            << "\n";
  std::cout << "DAE:       acc=" << dae.accepted << " rej=" << dae.rejected << " us=" << dae_us << " err=" << dae_err
            << "\n";
  std::cout << "DAE+Schur: acc=" << dae_schur.accepted << " rej=" << dae_schur.rejected << " us=" << schur_us
            << " err=" << schur_err << "\n";
  // Extended-family experiment: does constraining the peroxy radicals (the
  // next-fastest timescales in TS1's VOC chemistry) recover a step advantage?
  auto run_variant = [&](const std::string& label, const std::vector<std::string>& names)
  {
    std::vector<int> family;
    for (const auto& name : names)
      family.push_back(species_index.at(name));
    std::vector<double> family_initial = initial;
    Ts1QssaModel family_prototype(table, mech.species, family, std::vector<double>(family.size(), 1.0));
    family_prototype.Project(family_initial);
    std::vector<double> family_scale;
    for (std::size_t c = 0; c < family.size(); ++c)
    {
      double loss = 0.0;
      for (const std::size_t j : family_prototype.relevant_[c])
      {
        const auto& r = family_prototype.reactions_[j];
        const double stoich = NetStoichD(r, family[c]);
        if (stoich >= 0.0)
          continue;
        double rate = r.k;
        for (int q : r.reactants)
          rate *= family_initial[q];
        loss += -stoich * rate;
      }
      family_scale.push_back(1.0 / std::max(loss, 1e-300));
    }
    Ts1QssaModel family_constraint(table, mech.species, family, family_scale);
    auto family_solver = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>(options)
                             .SetSystem(micm::System(gas))
                             .SetReactions(processes)
                             .AddExternalModel(family_constraint)
                             .SetReorderState(false)
                             .Build();
    auto run = RunCase(family_solver, mech, family_initial, rtol, output_times);
    if (!run.converged)
    {
      std::cout << label << " (" << names.size() << " radicals): did not converge\n";
      return;
    }
    std::vector<double> samples;
    for (int rep = 0; rep < reps; ++rep)
      samples.push_back(TimedRun(family_solver, mech, family_initial, rtol, output_times));
    const double us = median(samples);
    const double err = MaxRelError(run, reference, output_times, t_skip);
    csv << label << ',' << rtol << ',' << run.accepted << ',' << run.rejected << ',' << us << ',' << err << '\n';
    std::cout << label << " (" << names.size() << " radicals): acc=" << run.accepted << " rej=" << run.rejected
              << " us=" << us << " err=" << err << "\n";
  };
  run_variant("dae_qssa9", { "O1D", "O", "OH", "HO2", "CH3O2", "C2H5O2", "CH3CO3", "C3H7O2", "PO2" });
  run_variant(
      "dae_qssa17",
      { "O1D",    "O",     "OH",    "HO2",   "CH3O2", "C2H5O2", "CH3CO3", "C3H7O2", "PO2",
        "MCO3",   "EO2",   "ALKO2", "MEKO2", "ENEO2", "XO2",    "RO2",    "HOCH2OO" });
  std::cout << "wrote ts1_dae.csv\n";
  return 0;
}
