#include "run_solver.hpp"

#include <micm/jit/jit_compiler.hpp>
#include <micm/jit/solver/jit_rosenbrock.hpp>
#include <micm/jit/solver/jit_solver_builder.hpp>
#include <micm/jit/solver/jit_solver_parameters.hpp>
#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/rate_constant/user_defined_rate_constant.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <gtest/gtest.h>
#include <omp.h>

using namespace micm;

template<std::size_t L>
using JitBuilder = JitSolverBuilder<JitRosenbrockSolverParameters, L>;

TEST(OpenMP, JITOneSolverManyStates)
{
  constexpr size_t n_threads = 8;

  auto a = micm::Species("A");
  auto b = micm::Species("B");
  auto c = micm::Species("C");

  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ a, b, c } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
                         .SetReactants({ a })
                         .SetProducts({ Yield(b, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r1" }))
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r2 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b, b })
                         .SetProducts({ Yield(b, 1), Yield(c, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r2" }))
                         .SetPhase(gas_phase)
                         .Build();

  micm::Process r3 = micm::ChemicalReactionBuilder()
                         .SetReactants({ b, c })
                         .SetProducts({ Yield(a, 1), Yield(c, 1) })
                         .SetRateConstant(micm::UserDefinedRateConstant({ .label_ = "r3" }))
                         .SetPhase(gas_phase)
                         .Build();

  auto reactions = std::vector<micm::Process>{ r1, r2, r3 };
  auto chemical_system = micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase });

  std::vector<std::vector<double>> results(n_threads);

  auto solver = JitBuilder<1>(RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(chemical_system)
                    .SetReactions(reactions)
                    .Build();

#pragma omp parallel num_threads(n_threads)
  {
    auto state = solver.GetState();
    std::vector<double> result = run_solver_on_thread_with_own_state(solver, state);
    results[omp_get_thread_num()] = result;
#pragma omp barrier
  }

  // compare each thread to thread 1
  for (int i = 1; i < n_threads; ++i)
  {
    for (int j = 0; j < results[0].size(); ++j)
    {
      EXPECT_EQ(results[0][j], results[i][j]);
    }
  }
}
