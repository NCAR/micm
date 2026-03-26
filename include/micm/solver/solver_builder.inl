// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
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
  template<class ExternalModel>
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
      StatePolicy>::AddExternalModelProcesses(ExternalModel&& model)
  {
    external_models_.emplace_back(
        ExternalModelProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>{ std::forward<decltype(model)>(model) });
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
      StatePolicy>::SetConstraints(std::vector<Constraint>&& constraints)
  {
    constraints_ = std::move(constraints);
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
      StatePolicy>::UnusedSpeciesCheck(const RatesPolicy& rates) const
  {
    if (ignore_unused_species_)
    {
      return;
    }

    auto used_species = rates.SpeciesUsed(reactions_);
    // Include species referenced by constraints (dependencies and algebraic targets)
    for (const auto& constraint : constraints_)
    {
      for (const auto& dep : constraint.SpeciesDependencies())
        used_species.insert(dep);
      used_species.insert(constraint.AlgebraicSpecies());
    }

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
      throw MicmException(MicmSeverity::Warning, MICM_ERROR_CATEGORY_SOLVER, MICM_SOLVER_ERROR_CODE_UNUSED_SPECIES, err_msg);
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
  inline std::unordered_map<std::string, std::size_t> SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>::GetSpeciesMap() const
  {
    std::unordered_map<std::string, std::size_t> species_map;
    std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> state_reordering;
    std::size_t index = 0;
    for (auto& name : system_.UniqueNames())
      species_map[name] = index++;

    if (reorder_state_)
    {
      // get unsorted Jacobian non-zero elements
      auto unsorted_rates = RatesPolicy(reactions_, species_map, external_models_);
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
  inline std::unordered_map<std::string, std::size_t> SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy>::GetCustomParameterMap() const
  {
    std::unordered_map<std::string, std::size_t> params{};
    std::vector<std::string> duplicates;

    auto add_param = [&params](const std::string& label, const std::string& source)
    {
      auto [it, added] = params.emplace_back(label, params.size());
      if (!added)
        duplicates.push_back(label + " (from " + source + ")");
    };

    // Include custom parameter labels from chemical reactions
    for (const auto& reaction : reactions_)
    {
      if (auto* process = std::get_if<ChemicalReaction>(&reaction.process_))
      {
        for (auto& label : process->rate_constant_->CustomParameters())
          add_params(label, "reaction");
      }
    }

    // Include custom parameter labels from external models
    for (const auto& model : system_.external_models_)
    {
      auto param_names = model.parameter_names_func_();
      for (const auto& label : param_names)
        add_param(label, "external_model");
    }

    // Include custom parameter labels from constraints
    for (const auto& constraint : constraints_)
    {
      auto param_names = constraint.GetParameterNames();
      for (const auto& label : param_names)
        add_param(label, "constraint");
    }

    if (!duplicates.empty())
    {
      std::ostringstream oss;
      oss << "Duplicate parameter labels detected:\n";
      for (const auto& d : duplicates)
          oss << "  - " << d << "\n";

      throw std::invalid_argument(oss.str());
    }

    return params;
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
      SetAbsoluteTolerances(std::vector<double>& tolerances, const std::unordered_map<std::string, std::size_t>& species_map)
          const
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
    if (!valid_system_)
    {
      throw MicmException(
          MicmSeverity::Error,
          MICM_ERROR_CATEGORY_SOLVER,
          MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SYSTEM,
          "Missing chemical system.");
    }

    std::size_t number_of_species = this->system_.StateSize();
    if (number_of_species == 0)
    {
      throw MicmException(
          MicmSeverity::Error,
          MICM_ERROR_CATEGORY_SOLVER,
          MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SPECIES,
          "Provided chemical system contains no species.");
    }

    using ConstraintSetPolicy = ConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>;
    using SolverPolicy =
        typename SolverParametersPolicy::template SolverType<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy>;

    auto species_map = this->GetSpeciesMap();

    RatesPolicy rates(reactions_, species_map, external_models_);

    this->UnusedSpeciesCheck(rates);
    auto nonzero_elements = rates.NonZeroJacobianElements();

    ConstraintSetPolicy constraint_set;
    std::set<std::size_t> algebraic_variable_ids;

    if (!constraints_.empty())
    {
      // Constraints replace selected species rows in the mass-matrix DAE formulation.
      // Pass species_map so constraints can resolve dependencies.
      constraint_set = ConstraintSetPolicy(constraints_, species_map);
      // Must set unqiue parameter names before the builder creates the parameter map.
      constraint_set.SetUniqueParameterNames();
      algebraic_variable_ids = constraint_set.AlgebraicVariableIds();
      rates.SetAlgebraicVariableIds(algebraic_variable_ids);

      // Filter kinetic sparsity entries from algebraic rows (they will be entirely replaced by constraints)
      for (auto it = nonzero_elements.begin(); it != nonzero_elements.end();)
      {
        if (algebraic_variable_ids.count(it->first) > 0)
          it = nonzero_elements.erase(it);
        else
          ++it;
      }

      // Merge constraint Jacobian elements with ODE Jacobian elements
      auto constraint_jac_elements = constraint_set.NonZeroJacobianElements();
      nonzero_elements.insert(constraint_jac_elements.begin(), constraint_jac_elements.end());
    }


    // The actual number of grid cells is not needed to construct the various solver objects
    auto jacobian = BuildJacobian<SparseMatrixPolicy>(nonzero_elements, 1, number_of_species, true);

    LinearSolverPolicy linear_solver(jacobian, 0);
    if constexpr (LuDecompositionInPlaceConcept<LuDecompositionPolicy, SparseMatrixPolicy>)
    {
      auto lu = LuDecompositionPolicy::template GetLUMatrix<SparseMatrixPolicy>(jacobian, 0, true);
      jacobian = std::move(lu);
    }

    // constraint set must ensure parameter names are unique before creating params map
    auto params_map = this->GetCustomParameterMap();

    rates.SetJacobianFlatIds(jacobian);
    rates.SetExternalModelFunctions(params_map, species_map, jacobian);

    if (constraint_set.Size() > 0)
    {
      constraint_set.SetJacobianFlatIds(jacobian);
      // TODO
      constraint_set.SetConstraintFunctions(params_map, species_map, jacobian);
    }

    std::vector<std::string> variable_names{ number_of_species };
    for (auto& species_pair : species_map)
      variable_names[species_pair.second] = species_pair.first;
    std::vector<std::string> labels{ params_map.size() };
    for (auto& param_pair : params_map)
      labels[param_pair.second] = param_pair.first;

    // Build mass-matrix diagonal: species rows default to ODE (1), rows replaced by constraints are algebraic (0).
    std::vector<double> mass_matrix_diagonal(number_of_species, 1.0);
    for (const auto variable_id : algebraic_variable_ids)
    {
      mass_matrix_diagonal[variable_id] = 0.0;
    }

    StateParameters state_parameters = { .number_of_species_ = number_of_species,
                                         .number_of_constraints_ = constraints_.size(),
                                         .number_of_rate_constants_ = this->reactions_.size(),
                                         .variable_names_ = variable_names,
                                         .custom_rate_parameter_labels_ = labels,
                                         .nonzero_jacobian_elements_ = nonzero_elements,
                                         .mass_matrix_diagonal_ = mass_matrix_diagonal };

    this->SetAbsoluteTolerances(state_parameters.absolute_tolerance_, species_map);

    // Create vector of functions to update external model state parameters
    std::vector<std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>> update_funcs;
    for (const auto& model : external_models_)
    {
      update_funcs.push_back(model.update_state_parameters_function_(params_map));
    }

    // make a copy of the options so that the builder can be used repeatedly
    // this matters because the absolute tolerances must be set to match the system size, and that may change
    auto options = this->options_;

    return Solver<SolverPolicy, StatePolicy>(
        SolverPolicy(std::move(linear_solver), std::move(rates), std::move(constraint_set)),
        state_parameters,
        options,
        reactions_,
        system_,
        update_funcs);
  }

}  // namespace micm
