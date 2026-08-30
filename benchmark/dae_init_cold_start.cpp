// Cold-start behaviour of DAE constraint initialization: does the damped Newton
// line search widen the basin of the consistent-initial-condition projection?
//
// Mechanism is the musica#956 Case-2 knife edge, unchanged from
// benchmark/musica956_knife_edge_main_api.cpp:
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
    constraints.emplace_back(EquilibriumConstraint<DenseMatrix, SparseMatrixT>(
        "henry",
        kAq,
        std::vector<StoichSpecies>{ { kGas, 1.0 } },
        std::vector<StoichSpecies>{ { kAq, 1.0 } },
        { .K_HLC_ref_ = 0.5, .delta_H_ = -60000.0 }));
    constraints.emplace_back(EquilibriumConstraint<DenseMatrix, SparseMatrixT>(
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
    double micros_;
  };

  Outcome Project(micm::Index max_backtracks, double temperature, micm::Real initial_xm, int repeats = 1)
  {
    auto params = RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    params.constraint_init_max_backtracks_ = max_backtracks;

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
      state.variables_[0][state.variable_map_.at("AQ")] = 1.0;
      state.variables_[0][state.variable_map_.at("HP")] = 1.0e3;
      state.variables_[0][state.variable_map_.at("XM")] = initial_xm;
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
    }

    std::sort(samples.begin(), samples.end());
    out.micros_ = samples[samples.size() / 2];
    return out;
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
      for (int i_temperature = 0; i_temperature <= 14; ++i_temperature)
      {
        const double T = 278.0 + static_cast<double>(i_temperature);
        for (const auto xm0 : kInitialXm)
        {
          const auto o = Project(variant.backtracks_, T, xm0);
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
    csv << "variant,start,temperature,initial_xm,status,constraint_init_iterations,jacobian_updates,"
           "decompositions,solves,median_us\n";
    csv << std::scientific << std::setprecision(17);

    // Warm: inside the full-step basin already. Consistent: on the manifold to
    // machine precision, so the projection should short-circuit before factoring.
    const std::vector<std::pair<const char*, micm::Real>> starts{ { "warm", 900.0 }, { "consistent", 0.0 } };

    for (const auto& variant : kVariants)
    {
      for (const auto& [label, xm0] : starts)
      {
        micm::Real start_xm = xm0;
        if (std::string(label) == "consistent")
        {
          // Converge first, then re-project from the converged value: the second
          // projection is the already-consistent common path R7 talks about.
          const auto converged = Project(24, 286.0, 900.0);
          start_xm = converged.xm_;
        }
        const auto o = Project(variant.backtracks_, 286.0, start_xm, 200);
        csv << variant.name_ << ',' << label << ',' << 286.0 << ',' << start_xm << ','
            << SolverStateToString(o.status_) << ',' << o.iterations_ << ',' << o.jacobians_ << ','
            << o.decompositions_ << ',' << o.solves_ << ',' << o.micros_ << '\n';
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
