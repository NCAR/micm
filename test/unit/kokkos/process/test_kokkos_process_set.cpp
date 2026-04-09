#include <gtest/gtest.h>
#include <micm/kokkos/process/kokkos_process_set.hpp>
#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/process/process.hpp>
#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>

TEST(KokkosProcessSet, ForcingTerms)
{
  micm::Species a("A"), b("B"), c("C");
  micm::Phase gas_phase{ "gas", std::vector<micm::PhaseSpecies>{ a, b, c } };

  micm::Process r1 = micm::ChemicalReactionBuilder()
    .SetReactants({ a, b })
    .SetProducts({ micm::StoichSpecies(c, 1.0) })
    .SetRateConstant(micm::ArrheniusRateConstant({ .A_ = 2.0 }))
    .SetPhase(gas_phase)
    .Build();

  std::vector<micm::Process> processes{ r1 };
  std::unordered_map<std::string, std::size_t> variable_map{ {"A", 0}, {"B", 1}, {"C", 2} };

  micm::KokkosProcessSet<micm::KokkosDenseMatrix<double>, micm::KokkosSparseMatrix<double>> process_set(processes, variable_map);

  micm::KokkosDenseMatrix<double> state_variables(1, 3);
  state_variables[0][0] = 1.0;
  state_variables[0][1] = 0.5;
  state_variables[0][2] = 0.0;
  state_variables.CopyToDevice();

  micm::KokkosDenseMatrix<double> rate_constants(1, 1);
  rate_constants[0][0] = 2.0;
  rate_constants.CopyToDevice();

  micm::KokkosDenseMatrix<double> forcing(1, 3, 0.0);
  forcing.CopyToDevice();

  // Make sure to sync all data to device
  state_variables.CopyToDevice();
  rate_constants.CopyToDevice();
  forcing.CopyToDevice();
  Kokkos::fence();

  micm::State<micm::KokkosDenseMatrix<double>, micm::KokkosSparseMatrix<double>> state(
    micm::StateParameters{ .number_of_rate_constants_ = 1, .variable_names_{ "A", "B", "C" } }, 1);
  state.rate_constants_ = rate_constants;
  state.rate_constants_.CopyToDevice();

  process_set.AddForcingTerms(state, state_variables, forcing);

  forcing.CopyToHost();


  // Rate = 2.0 * 1.0 * 0.5 = 1.0
  // d[A]/dt = -1.0
  // d[B]/dt = -1.0
  // d[C]/dt = 1.0
  EXPECT_NEAR(forcing[0][0], -1.0, 1e-8);
  EXPECT_NEAR(forcing[0][1], -1.0, 1e-8);
  EXPECT_NEAR(forcing[0][2], 1.0, 1e-8);
}
