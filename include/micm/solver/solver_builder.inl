// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

enum class MicmSolverBuilderErrc
{
  UnusedSpecies = 1,  // Unused species present in the chemical system
  MissingChemicalSystem = 2, // Missing chemical system
  MissingReactions = 3, // Missing processes
  MissingChemicalSpecies = 4 // Missing chemical species
};

namespace std
{
  template<>
  struct is_error_condition_enum<MicmSolverBuilderErrc> : true_type
  {
  };
}  // namespace std

namespace
{
  class SolverBuilderErrorCategory : public std::error_category
  {
   public:
    const char* name() const noexcept override
    {
      return "MICM Solver Builder";
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmSolverBuilderErrc>(ev))
      {
        case MicmSolverBuilderErrc::UnusedSpecies:
          return "Unused species present in the chemical system. Use the ignore_unused_species_ parameter to allow unused "
                 "species in the solve.";
        case MicmSolverBuilderErrc::MissingChemicalSystem:
          return "Missing chemical system. Use the SetSystem function to set the chemical system.";
        case MicmSolverBuilderErrc::MissingReactions:
          return "Missing reactions. Use the SetReactions function to set the processes.";
        case MicmSolverBuilderErrc::MissingChemicalSpecies:
          return "Provided chemical system contains no species.";
        default: return "Unknown error";
      }
    }
  };

  const SolverBuilderErrorCategory solverBuilderErrorCategory{};
}  // namespace

inline std::error_code make_error_code(MicmSolverBuilderErrc e)
{
  return { static_cast<int>(e), solverBuilderErrorCategory };
}

namespace micm
{
  template<class SolverParametersPolicy, class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy, class LinearSolverPolicy>
  inline SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>&
  SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>::SetSystem(const System& system)
  {
    system_ = system;
    valid_system_ = true;
    return *this;
  }

  template<class SolverParametersPolicy, class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy, class LinearSolverPolicy>
  inline SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>&
  SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>::SetReactions(const std::vector<Process>& reactions)
  {
    reactions_ = reactions;
    valid_reactions_ = reactions_.size() > 0;
    return *this;
  }

  template<class SolverParametersPolicy, class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy, class LinearSolverPolicy>
  inline SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>&
  SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>::SetNumberOfGridCells(int number_of_grid_cells)
  {
    number_of_grid_cells_ = number_of_grid_cells;
    return *this;
  }

  template<class SolverParametersPolicy, class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy, class LinearSolverPolicy>
  inline SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>&
  SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>::SetIgnoreUnusedSpecies(bool ignore_unused_species)
  {
    ignore_unused_species_ = ignore_unused_species;
    return *this;
  }

  template<class SolverParametersPolicy, class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy, class LinearSolverPolicy>
  inline SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>&
  SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>::SetReorderState(bool reorder_state)
  {
    reorder_state_ = reorder_state;
    return *this;
  }

  template<class SolverParametersPolicy, class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy, class LinearSolverPolicy>
  inline void SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>::UnusedSpeciesCheck()
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
      throw std::system_error(make_error_code(MicmSolverBuilderErrc::UnusedSpecies), err_msg);
    }
  }

  template<class SolverParametersPolicy, class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy, class LinearSolverPolicy>
  inline std::map<std::string, std::size_t> SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>::GetSpeciesMap() const
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

  template<class SolverParametersPolicy, class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy, class LinearSolverPolicy>
  inline void SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>::SetAbsoluteTolerances(
      std::vector<double>& tolerances,
      const std::map<std::string, std::size_t>& species_map) const
  {
    // if the tolerances aren't already set, initialize them and then set based off of information in the system
    if (tolerances.size() != species_map.size())
    {
      tolerances = std::vector<double>(species_map.size(), 1e-3);
      for (auto& species : system_.gas_phase_.species_)
      {
        if (species.HasProperty("absolute tolerance"))
        {
          tolerances[species_map.at(species.name_)] = species.template GetProperty<double>("absolute tolerance");
        }
      }
      for (auto& phase : system_.phases_)
      {
        for (auto& species : phase.second.species_)
        {
          if (species.HasProperty("absolute tolerance"))
          {
            tolerances[species_map.at(species.name_)] = species.template GetProperty<double>("absolute tolerance");
          }
        }
      }
    }
  }

  template<class SolverParametersPolicy, class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy, class LinearSolverPolicy>
  inline std::vector<std::string> SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>::GetCustomParameterLabels() const
  {
    std::vector<std::string> param_labels{};
    for (const auto& reaction : reactions_)
      if (reaction.rate_constant_)
        for (auto& label : reaction.rate_constant_->CustomParameters())
          param_labels.push_back(label);
    return param_labels;
  }

  template<class SolverParametersPolicy, class DenseMatrixPolicy, class SparseMatrixPolicy, class RatesPolicy, class LinearSolverPolicy>
  inline auto SolverBuilder<SolverParametersPolicy, DenseMatrixPolicy, SparseMatrixPolicy, RatesPolicy, LinearSolverPolicy>::Build()
  {
    if (!valid_system_)
    {
      throw std::system_error(make_error_code(MicmSolverBuilderErrc::MissingChemicalSystem), "Missing chemical system.");
    }
    if (!valid_reactions_)
    {
      throw std::system_error(make_error_code(MicmSolverBuilderErrc::MissingReactions), "Missing reactions.");
    }
    using SolverPolicy = typename SolverParametersPolicy::template SolverType<RatesPolicy, LinearSolverPolicy>;
    auto species_map = this->GetSpeciesMap();
    auto labels = this->GetCustomParameterLabels();
    std::size_t number_of_species = this->system_.StateSize();
    if (number_of_species == 0)
    {
      throw std::system_error(make_error_code(MicmSolverBuilderErrc::MissingChemicalSpecies), "Provided chemical system contains no species.");
    }

    this->UnusedSpeciesCheck();
    this->SetAbsoluteTolerances(this->options_.absolute_tolerance_, species_map);

    RatesPolicy rates(this->reactions_, species_map);
    auto nonzero_elements = rates.NonZeroJacobianElements();
    auto jacobian = BuildJacobian<SparseMatrixPolicy>(nonzero_elements, this->number_of_grid_cells_, number_of_species);

    rates.SetJacobianFlatIds(jacobian);
    LinearSolverPolicy linear_solver(jacobian, 1e-30);

    std::vector<std::string> variable_names{ number_of_species };
    for (auto& species_pair : species_map)
      variable_names[species_pair.second] = species_pair.first;

    StateParameters state_parameters = { .number_of_grid_cells_ = this->number_of_grid_cells_,
                                          .number_of_species_ = number_of_species,
                                          .number_of_rate_constants_ = this->reactions_.size(),
                                          .variable_names_ = variable_names,
                                          .custom_rate_parameter_labels_ = labels,
                                          .nonzero_jacobian_elements_ = nonzero_elements };

    return Solver<SolverPolicy, State<DenseMatrixPolicy, SparseMatrixPolicy>>(
          SolverPolicy(
              this->options_, std::move(linear_solver), std::move(rates), jacobian),
          state_parameters,
          this->number_of_grid_cells_,
          number_of_species,
          this->reactions_.size(),
          this->reactions_);
  }

}  // namespace micm