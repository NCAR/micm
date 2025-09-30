// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace micm
{
  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy>
  inline SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>&
  SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>::SetSystem(const System& system)
  {
    system_ = system;
    valid_system_ = true;
    return *this;
  }

  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy>
  inline SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>&
  SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>::SetReactions(const std::vector<Process>& reactions)
  {
    reactions_ = reactions;
    valid_reactions_ = reactions_.size() > 0;
    return *this;
  }

  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy>
  inline SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>&
  SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>::SetIgnoreUnusedSpecies(bool ignore_unused_species)
  {
    ignore_unused_species_ = ignore_unused_species;
    return *this;
  }

  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy>
  inline SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>&
  SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>::SetReorderState(bool reorder_state)
  {
    reorder_state_ = reorder_state;
    return *this;
  }

  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy>
  inline void SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>::UnusedSpeciesCheck() const
  {
    if (ignore_unused_species_)
    {
      return;
    }

    auto used_species = RatesPolicy::SpeciesUsed(reactions_);
    auto available_species = system_.UniqueNames();
    std::sort(available_species.begin(), available_species.end());
    std::set<std::string> unused_species;
    std::set_difference(
        available_species.begin(),
        available_species.end(),
        used_species.begin(),
        used_species.end(),
        std::inserter(unused_species, unused_species.begin()));
    if (unused_species.size() > 0)
    {
      std::string err_msg = "Unused species in chemical system:";
      for (auto& species : unused_species)
        err_msg += " '" + species + "'";
      err_msg += ".";
      throw std::system_error(make_error_code(MicmSolverErrc::UnusedSpecies), err_msg);
    }
  }

  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy>
  inline std::map<std::string, std::size_t> SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>::GetSpeciesMap() const
  {
    std::map<std::string, std::size_t> species_map;
    std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> state_reordering;
    std::size_t index = 0;
    for (auto& name : system_.UniqueNames())
      species_map[name] = index++;

    if (reorder_state_)
    {
      // get unsorted Jacobian non-zero elements
      auto unsorted_rates = RatesPolicy(reactions_, species_map);
      auto unsorted_jac_elements = unsorted_rates.NonZeroJacobianElements();

      using Matrix = typename DenseMatrixPolicy::IntMatrix;
      Matrix unsorted_jac_non_zeros(system_.StateSize(), system_.StateSize(), 0);
      for (auto& elem : unsorted_jac_elements)
        unsorted_jac_non_zeros[elem.first][elem.second] = 1;
      auto reorder_map = DiagonalMarkowitzReorder<Matrix>(unsorted_jac_non_zeros);

      state_reordering = [=](const std::vector<std::string>& variables, const std::size_t i)
      { return variables[reorder_map[i]]; };

      index = 0;
      for (auto& name : system_.UniqueNames(state_reordering))
        species_map[name] = index++;
    }

    return species_map;
  }

  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy>
  inline void SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>::
      SetAbsoluteTolerances(std::vector<double>& tolerances, const std::map<std::string, std::size_t>& species_map) const
  {
    tolerances = std::vector<double>(species_map.size(), 1e-3);
    for (auto& phase_species : system_.gas_phase_.phase_species_)
    {
      auto& species = phase_species.species_;
      if (species.HasProperty("absolute tolerance"))
      {
        tolerances[species_map.at(species.name_)] = species.template GetProperty<double>("absolute tolerance");
      }
    }
    for (auto& phase : system_.phases_)
    {
      for (auto& phase_species : phase.second.phase_species_)
      {
        auto& species = phase_species.species_;
        if (species.HasProperty("absolute tolerance"))
        {
          tolerances[species_map.at(species.name_)] = species.template GetProperty<double>("absolute tolerance");
        }
      }
    }
  }

  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy>
  inline std::vector<std::string> SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>::GetCustomParameterLabels() const
  {
    std::vector<std::string> param_labels{};

    for (const auto& reaction : reactions_)
    {
      if (auto* process = std::get_if<ChemicalReaction>(&reaction.process_))
      {
        for (auto& label : process->rate_constant_->CustomParameters())
        {
          param_labels.push_back(label);
        }
      }
    }
    return param_labels;
  }

  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy>
  inline auto SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>::Build() const
  {
    // make a copy of the options so that the builder can be used repeatedly
    // this matters because the absolute tolerances must be set to match the system size, and that may change
    auto options = this->options_;

    if (!valid_system_)
    {
      throw std::system_error(make_error_code(MicmSolverErrc::MissingChemicalSystem), "Missing chemical system.");
    }
    if (!valid_reactions_)
    {
      throw std::system_error(make_error_code(MicmSolverErrc::MissingProcesses), "Missing processes.");
    }
    using SolverPolicy = typename SolverParametersPolicy::template SolverType<RatesPolicy, LinearSolverPolicy>;
    auto species_map = this->GetSpeciesMap();
    auto labels = this->GetCustomParameterLabels();
    std::size_t number_of_species = this->system_.StateSize();
    if (number_of_species == 0)
    {
      throw std::system_error(
          make_error_code(MicmSolverErrc::MissingChemicalSpecies), "Provided chemical system contains no species.");
    }

    this->UnusedSpeciesCheck();

    RatesPolicy rates(this->reactions_, species_map);
    auto nonzero_elements = rates.NonZeroJacobianElements();
    // The actual number of grid cells is not needed to construct the various solver objects
    auto jacobian = BuildJacobian<SparseMatrixPolicy>(nonzero_elements, 1, number_of_species, true);

    LinearSolverPolicy linear_solver(jacobian, 0);
    if constexpr (LuDecompositionInPlaceConcept<LuDecompositionPolicy, SparseMatrixPolicy>)
    {
      auto lu = LuDecompositionPolicy::template GetLUMatrix<SparseMatrixPolicy>(jacobian, 0, true);
      jacobian = std::move(lu);
    }
    rates.SetJacobianFlatIds(jacobian);

    std::vector<std::string> variable_names{ number_of_species };
    for (auto& species_pair : species_map)
      variable_names[species_pair.second] = species_pair.first;

    StateParameters state_parameters = { .number_of_species_ = number_of_species,
                                         .number_of_rate_constants_ = this->reactions_.size(),
                                         .variable_names_ = variable_names,
                                         .custom_rate_parameter_labels_ = labels,
                                         .nonzero_jacobian_elements_ = nonzero_elements };

    this->SetAbsoluteTolerances(state_parameters.absolute_tolerance_, species_map);

    return Solver<SolverPolicy, StatePolicy>(
        SolverPolicy(std::move(linear_solver), std::move(rates), jacobian, number_of_species),
        state_parameters,
        options,
        this->reactions_,
        this->system_);
  }

}  // namespace micm