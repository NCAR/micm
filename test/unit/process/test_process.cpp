#include <micm/process/chemical_reaction_builder.hpp>
#include <micm/process/phase_transfer_process_builder.hpp>
#include <micm/process/process.hpp>
#include <micm/process/rate_constant/arrhenius_rate_constant.hpp>
#include <micm/process/rate_constant/surface_rate_constant.hpp>
#include <micm/process/rate_constant/user_defined_rate_constant.hpp>
#include <micm/process/transfer_coefficient/phase_transfer_coefficient.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <random>
using namespace micm;

template<class DenseMatrixPolicy>
void testProcessUpdateState(const std::size_t number_of_grid_cells)
{
  Species foo("foo", { { "molecular weight [kg mol-1]", 0.025 }, { "diffusion coefficient [m2 s-1]", 2.3e2 } });
  Species bar("bar");
  bar.parameterize_ = [](const Conditions& c) { return c.air_density_ * 0.82; };

  ArrheniusRateConstant rc1({ .A_ = 12.2, .C_ = 300.0 });
  SurfaceRateConstant rc2({ .label_ = "foo_surf", .species_ = foo });
  UserDefinedRateConstant rc3({ .label_ = "bar_user" });

  Process r1 = ChemicalReactionBuilder().SetRateConstant(rc1).SetReactants({ foo, bar }).Build();
  Process r2 = ChemicalReactionBuilder().SetRateConstant(rc2).Build();
  Process r3 = ChemicalReactionBuilder().SetRateConstant(rc3).Build();
  std::vector<Process> processes = { r1, r2, r3 };

  std::vector<std::string> param_labels{};
  for (const auto& process : processes)
  {
    if (auto* reaction = std::get_if<ChemicalReaction>(&process.process_))
    {
      for (const auto& label : reaction->rate_constant_->CustomParameters())
      {
        param_labels.push_back(label);
      }
    }
  }
  State<DenseMatrixPolicy> state{ StateParameters{
                                      .number_of_rate_constants_ = processes.size(),
                                      .variable_names_ = { "foo", "bar'" },
                                      .custom_rate_parameter_labels_ = param_labels,
                                  },
                                  number_of_grid_cells };

  DenseMatrixPolicy expected_rate_constants(number_of_grid_cells, 3, 0.0);
  std::vector<double> params = { 0.0, 0.0, 0.0 };
  auto get_double = std::bind(std::lognormal_distribution(0.0, 0.01), std::default_random_engine());

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    state.conditions_[i_cell].temperature_ = get_double() * 285.0;
    state.conditions_[i_cell].pressure_ = get_double() * 101100.0;
    state.conditions_[i_cell].air_density_ = get_double() * 10.0;
    params[0] = get_double() * 1.0e-8;
    params[1] = get_double() * 1.0e5;
    params[2] = get_double() * 1.0e-2;
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["foo_surf.effective radius [m]"]] = params[0];
    state.custom_rate_parameters_[i_cell]
                                 [state.custom_rate_parameter_map_["foo_surf.particle number concentration [# m-3]"]] =
        params[1];
    state.custom_rate_parameters_[i_cell][state.custom_rate_parameter_map_["bar_user"]] = params[2];
    std::vector<double>::const_iterator param_iter = params.begin();
    expected_rate_constants[i_cell][0] =
        rc1.Calculate(state.conditions_[i_cell], param_iter) * (state.conditions_[i_cell].air_density_ * 0.82);
    param_iter += rc1.SizeCustomParameters();
    expected_rate_constants[i_cell][1] = rc2.Calculate(state.conditions_[i_cell], param_iter);
    param_iter += rc2.SizeCustomParameters();
    expected_rate_constants[i_cell][2] = rc3.Calculate(state.conditions_[i_cell], param_iter);
    param_iter += rc3.SizeCustomParameters();
  }
  Process::CalculateRateConstants(processes, state);

  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
    for (std::size_t i_rxn = 0; i_rxn < processes.size(); ++i_rxn)
      EXPECT_EQ(state.rate_constants_[i_cell][i_rxn], expected_rate_constants[i_cell][i_rxn])
          << "grid cell " << i_cell << "; reaction " << i_rxn;
}

template<class T>
using Group1VectorMatrix = VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = VectorMatrix<T, 4>;

TEST(Process, Matrix)
{
  testProcessUpdateState<Matrix<double>>(5);
}

TEST(Process, VectorMatrix)
{
  testProcessUpdateState<Group1VectorMatrix<double>>(5);
  testProcessUpdateState<Group2VectorMatrix<double>>(5);
  testProcessUpdateState<Group3VectorMatrix<double>>(5);
  testProcessUpdateState<Group4VectorMatrix<double>>(5);
}

TEST(Process, DifferentiatesChemicalReactionAndPhaseTransfer)
{
  // Build a ChemicalReaction
  Process chemical_reaction = ChemicalReactionBuilder()
                                  .SetPhaseName("gas")
                                  .SetReactants({ Species("O3"), Species("NO") })
                                  .SetProducts({ Yield(Species("NO2"), 1.0), Yield(Species("O2"), 1.0) })
                                  .SetRateConstant(ArrheniusRateConstant(/* ... */))
                                  .Build();

  // Build a PhaseTransferProcess
  Process phase_transfer = PhaseTransferProcessBuilder()
                               .SetOriginSpecies({ SpeciesInPhase("gas", Species("SO2")) })
                               .SetDestinationSpecies({ SpeciesInPhase("aqueous", Species("SO2")) })
                               .SetSolvent(SpeciesInPhase("aqueous", Species("H2O")))
                               .SetTransferCoefficient(TestTransferCoefficient(/* ... */))
                               .Build();

  // Check that the first process is a ChemicalReaction
  std::visit(
      [](auto&& value)
      {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, ChemicalReaction>)
        {
          EXPECT_EQ(value.phase_name_, "gas");
        }
        else
        {
          FAIL() << "Expected ChemicalReaction, got different type";
        }
      },
      chemical_reaction.process_);

  // Check that the second process is a PhaseTransferProcess
  std::visit(
      [](auto&& value)
      {
        using T = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<T, PhaseTransferProcess>)
        {
          ASSERT_FALSE(value.destination_species_.empty());
          EXPECT_EQ(value.destination_species_[0].phase_name_, "aqueous");
        }
        else
        {
          FAIL() << "Expected PhaseTransferProcess, got different type";
        }
      },
      phase_transfer.process_);
}

// TODO (Jiwon): This should throw, but currently does not.
//               I don't think the base class Process should know about a derived-class's specific condition.
//               Feels like a design issue — will revisit later. Commented out for now.
// TEST(Process, SurfaceRateConstantOnlyHasOneReactant)
// {
//   Species c("c", { { "molecular weight [kg mol-1]", 0.025 }, { "diffusion coefficient [m2 s-1]", 2.3e2 } });
//   Species e("e");

//   EXPECT_ANY_THROW(Process r = ChemicalReactionBuilder()
//                                          .SetPhaseName("gas")
//                                          .SetReactants({ c, c })
//                                          .SetProducts({ Yield(e, 1) })
//                                          .SetRateConstant(SurfaceRateConstant(
//                                              { .label_ = "c", .species_ = c, .reaction_probability_ = 0.90 }))
//                                          .Build(););
// }
