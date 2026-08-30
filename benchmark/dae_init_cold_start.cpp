// Cold-start behaviour of DAE constraint initialization: does the damped Newton
// line search improve convergence reach within the configured update budget?
//
// Reduced mechanism inspired by the musica#956 Case-2 knife edge, unchanged from
// the specification's Henry/dissociation reproducer:
//   Henry's law   K1(T)*[GAS] = [AQ]      (Van't Hoff, exothermic)
//   Dissociation  K2*[AQ]     = [XM]^2    (quadratic in the solved ion)
//
// The projection is driven directly through InitializeConstraints rather than
// through Solve, so the measurement isolates the initialization and reports its
// work counters without the integration loop's contribution.
//
// Writes, in the current working directory:
//   dae_init_basin.csv -- D1: status over (initial XM) x temperature x variant
//   dae_init_cost.csv  -- D2: initialization work and wall time, warm/consistent
//
// Build: gated by MICM_ENABLE_DAE_BENCHMARKS (see benchmark/CMakeLists.txt).
#include <micm/CPU.hpp>
#include <micm/constraint/constraint.hpp>
#include <micm/constraint/types/equilibrium_constraint.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <ios>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace micm;

using StandardBuilder =
    CpuSolverBuilder<RosenbrockSolverParameters, Matrix<micm::Real>, SparseMatrix<micm::Real, SparseMatrixStandardOrdering>>;
using DenseMatrix = StandardBuilder::DenseMatrixPolicyType;
using SparseMatrixT = StandardBuilder::SparseMatrixPolicyType;

namespace
{
  const auto kGas = Species("GAS");
  const auto kAq = Species("AQ");
  const auto kHp = Species("HP");
  const auto kXm = Species("XM");
  const auto kSink = Species("SINK");

  Phase CloudPhase()
  {
    return Phase{ "cloud", std::vector<PhaseSpecies>{ kGas, kAq, kHp, kXm, kSink } };
  }

  Process LossReaction(const Phase& phase)
  {
    return ChemicalReactionBuilder()
        .SetReactants({ kGas })
        .SetProducts({ StoichSpecies(kSink, 1) })
        .SetRateConstant(ArrheniusRateConstantParameters{ .A_ = 1.0e-12 })
        .SetPhase(phase)
        .Build();
  }

  std::vector<Constraint<DenseMatrix, SparseMatrixT>> KnifeEdgeConstraints()
  {
    std::vector<Constraint<DenseMatrix, SparseMatrixT>> constraints;
    constraints.emplace_back(
        EquilibriumConstraint<DenseMatrix, SparseMatrixT>(
            "henry",
            kAq,
            std::vector<StoichSpecies>{ { kGas, 1.0 } },
            std::vector<StoichSpecies>{ { kAq, 1.0 } },
            { .K_HLC_ref_ = 0.5, .delta_H_ = -60000.0 }));
    constraints.emplace_back(
        EquilibriumConstraint<DenseMatrix, SparseMatrixT>(
            "dissoc",
            kXm,
            std::vector<StoichSpecies>{ { kAq, 1.0 } },
            std::vector<StoichSpecies>{ { kXm, 2.0 } },
            { .K_HLC_ref_ = 1.0, .delta_H_ = 0.0 }));
    return constraints;
  }

  /// @brief One projection: status, the post-init AQ and XM, and the work it cost.
  struct Outcome
  {
    SolverState status_;
    micm::Real aq_;
    micm::Real xm_;
    micm::Index iterations_;
    micm::Index jacobians_;
    micm::Index decompositions_;
    micm::Index solves_;
    micm::Index forcing_calls_;
    double micros_;
  };

  struct AlgebraicStart
  {
    micm::Real aq_;
    micm::Real xm_;
  };

  RosenbrockSolverParameters Parameters(micm::Index max_backtracks)
  {
    auto params = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    params.constraint_init_max_backtracks_ = max_backtracks;
    return params;
  }

  Outcome
  Project(const RosenbrockSolverParameters& params, double temperature, const AlgebraicStart& initial, int repeats = 1)
  {
    auto constraints = KnifeEdgeConstraints();
    const auto phase = CloudPhase();
    auto solver = StandardBuilder(params)
                      .SetSystem(System(phase))
                      .SetReactions({ LossReaction(phase) })
                      .SetConstraints(std::move(constraints))
                      .SetReorderState(false)
                      .Build();

    Outcome out{};
    std::vector<double> samples;
    samples.reserve(static_cast<std::size_t>(repeats));

    for (int r = 0; r < repeats; ++r)
    {
      auto state = solver.GetState(1);
      state.variables_[0][state.variable_map_.at("GAS")] = 4.0e6;
      state.variables_[0][state.variable_map_.at("AQ")] = initial.aq_;
      state.variables_[0][state.variable_map_.at("HP")] = 1.0e3;
      state.variables_[0][state.variable_map_.at("XM")] = initial.xm_;
      state.variables_[0][state.variable_map_.at("SINK")] = 0.0;
      state.conditions_[0].temperature_ = temperature;
      state.conditions_[0].pressure_ = 85000.0;
      solver.UpdateStateParameters(state);
      state.variables_.CopyToDevice();

      SolverStats stats;
      const auto t0 = std::chrono::steady_clock::now();
      const auto status = solver.solver_.InitializeConstraints(state, params, stats);
      const auto t1 = std::chrono::steady_clock::now();
      state.variables_.CopyToHost();

      samples.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
      out.status_ = status;
      out.aq_ = state.variables_[0][state.variable_map_.at("AQ")];
      out.xm_ = state.variables_[0][state.variable_map_.at("XM")];
      out.iterations_ = stats.constraint_init_iterations_;
      out.jacobians_ = stats.jacobian_updates_;
      out.decompositions_ = stats.decompositions_;
      out.solves_ = stats.solves_;
      // InitializeConstraints deliberately does not add constraint residual evaluations to
      // SolverStats::function_calls_: changing that counter would violate the exact #1083 opt-out.
      // In these finite benchmark paths, every Newton pass evaluates the current residual once and
      // every solve beyond the current-Jacobian solve corresponds to one candidate residual. The
      // exact-zero path has one current-residual evaluation and no solves.
      out.forcing_calls_ = stats.constraint_init_iterations_ + stats.solves_ - stats.jacobian_updates_;
    }

    std::sort(samples.begin(), samples.end());
    out.micros_ = samples[samples.size() / 2];
    return out;
  }

  AlgebraicStart ExactConsistentStart(double temperature)
  {
    auto exact_params = Parameters(24);
    exact_params.constraint_init_tolerance_ = 0.0;
    exact_params.constraint_init_max_iterations_ = 50;

    const auto converged = Project(exact_params, temperature, { 1.0, 900.0 });
    if (converged.status_ != SolverState::Converged)
    {
      throw std::runtime_error("could not construct an exactly consistent DAE state");
    }

    // Prove that rebuilding the same mechanism from the saved complete algebraic vector takes the
    // exact-zero-residual short circuit: one residual evaluation, no Jacobian, factorization, or solve.
    const auto verified = Project(Parameters(0), temperature, { converged.aq_, converged.xm_ });
    if (verified.status_ != SolverState::Converged || verified.iterations_ != 1 || verified.forcing_calls_ != 1 ||
        verified.jacobians_ != 0 || verified.decompositions_ != 0 || verified.solves_ != 0)
    {
      throw std::runtime_error("constructed DAE state does not take the exact-zero-residual short circuit");
    }

    return { converged.aq_, converged.xm_ };
  }

  const std::vector<micm::Real> kInitialXm{ 1.0, 10.0, 100.0, 300.0, 600.0, 800.0, 900.0, 1200.0, 2000.0 };

  struct Variant
  {
    const char* name_;
    micm::Index backtracks_;
  };
  const std::vector<Variant> kVariants{ { "undamped", 0 }, { "damped", 24 } };

  void RunBasinMap()
  {
    std::ofstream csv("dae_init_basin.csv");
    csv << "variant,temperature,initial_xm,status,final_aq,final_xm,constraint_init_iterations,solves\n";
    csv << std::scientific << std::setprecision(17);

    for (const auto& variant : kVariants)
    {
      const auto params = Parameters(variant.backtracks_);
      for (int i_temperature = 0; i_temperature <= 14; ++i_temperature)
      {
        const double T = 278.0 + static_cast<double>(i_temperature);
        for (const auto xm0 : kInitialXm)
        {
          const auto o = Project(params, T, { 1.0, xm0 });
          csv << variant.name_ << ',' << T << ',' << xm0 << ',' << SolverStateToString(o.status_) << ',' << o.aq_ << ','
              << o.xm_ << ',' << o.iterations_ << ',' << o.solves_ << '\n';
        }
      }
    }
    std::printf("wrote dae_init_basin.csv\n");
  }

  void RunCost()
  {
    std::ofstream csv("dae_init_cost.csv");
    csv << "variant,start,temperature,initial_aq,initial_xm,status,constraint_init_iterations,forcing_calls,"
           "jacobian_updates,decompositions,solves,median_us\n";
    csv << std::scientific << std::setprecision(17);

    // Warm: inside the full-step basin already. Consistent: the complete algebraic vector has an
    // exactly zero residual and was independently verified to short-circuit before factorization.
    const std::vector<std::pair<const char*, AlgebraicStart>> starts{ { "warm", { 1.0, 900.0 } },
                                                                      { "consistent", ExactConsistentStart(286.0) } };

    for (const auto& variant : kVariants)
    {
      const auto params = Parameters(variant.backtracks_);
      for (const auto& [label, initial] : starts)
      {
        const auto o = Project(params, 286.0, initial, 200);
        if (std::string(label) == "consistent" &&
            (o.iterations_ != 1 || o.forcing_calls_ != 1 || o.jacobians_ != 0 || o.decompositions_ != 0 || o.solves_ != 0))
        {
          throw std::runtime_error("consistent cost row performed work beyond its residual check");
        }
        csv << variant.name_ << ',' << label << ',' << 286.0 << ',' << initial.aq_ << ',' << initial.xm_ << ','
            << SolverStateToString(o.status_) << ',' << o.iterations_ << ',' << o.forcing_calls_ << ',' << o.jacobians_
            << ',' << o.decompositions_ << ',' << o.solves_ << ',' << o.micros_ << '\n';
      }
    }
    std::printf("wrote dae_init_cost.csv\n");
  }
}  // namespace

int main()
{
  RunBasinMap();
  RunCost();
  return 0;
}
