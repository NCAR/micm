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
    // Include species referenced by external model constraints
    for (const auto& model : external_constraint_models_)
    {
      auto deps = model.species_dependencies_func_();
      used_species.insert(deps.begin(), deps.end());
      auto alg_names = model.algebraic_variable_names_func_();
      used_species.insert(alg_names.begin(), alg_names.end());
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

    auto add_param = [&params, &duplicates](const std::string& label, const std::string& source)
    {
      auto [it, added] = params.emplace(label, params.size());
      if (!added)
        duplicates.push_back(label + " (from " + source + ")");
    };

    // Include custom parameter labels from chemical reactions
    for (const auto& reaction : reactions_)
    {
      if (auto* process = std::get_if<ChemicalReaction>(&reaction.process_))
      {
        for (auto& label : process->rate_constant_->CustomParameters())
          add_param(label, "reaction");
      }
    }

    // Include custom parameter labels from external models
    for (const auto& model : system_.external_models_)
    {
      auto param_names = model.parameter_names_func_();
      for (const auto& label : param_names)
        add_param(label, "external_model");
    }

    if (!duplicates.empty())
    {
      std::ostringstream oss;
      oss << "Duplicate parameter labels detected:\n";
      for (const auto& d : duplicates)
        oss << "  - " << d << "\n";

      throw MicmException(
          MicmSeverity::Error, MICM_ERROR_CATEGORY_SOLVER, MICM_SOLVER_ERROR_CODE_DUPLICATE_PARAMETER, oss.str());
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
      StatePolicy>::Build()
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

    // Constraints are not supported with CUDA matrix policies
    constexpr bool is_cuda_policy = requires(DenseMatrixPolicy m) { m.CopyToDevice(); m.CopyToHost(); };
    if constexpr (is_cuda_policy)
    {
      if (!constraints_.empty() || !external_constraint_models_.empty())
      {
        throw MicmException(
            MicmSeverity::Error,
            MICM_ERROR_CATEGORY_SOLVER,
            MICM_SOLVER_ERROR_CODE_CUDA_CONSTRAINTS_UNSUPPORTED,
            "Constraints are not supported with CUDA matrix policies.");
      }
    }

    using ConstraintSetPolicy = ConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>;
    using SolverPolicy =
        typename SolverParametersPolicy::template SolverType<RatesPolicy, LinearSolverPolicy, ConstraintSetPolicy>;

    // Build ProcessSet
    auto species_map = this->GetSpeciesMap();
    RatesPolicy rates(reactions_, species_map, external_models_);
    this->UnusedSpeciesCheck(rates);
    auto nonzero_elements = rates.NonZeroJacobianElements();

    auto params_map = this->GetCustomParameterMap();

    // Create vector of functions to update external model state parameters
    // (compiled after all params are added to params_map — see below)
    std::vector<std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>> update_state_param_funcs;

    // Build constraint set
    ConstraintSetPolicy constraint_set;

    // Build mass-matrix diagonal: species rows default to ODE (1), rows replaced by constraints are algebraic (0).
    std::vector<double> mass_matrix_diagonal(number_of_species, 1.0);
    std::set<std::size_t> algebraic_variable_ids;

    if (!constraints_.empty())
    {
      // Constraints replace selected species rows in the mass-matrix DAE formulation.
      // Pass species_map so constraints can resolve dependencies.
      constraint_set = ConstraintSetPolicy(std::move(constraints_), species_map);

      // Set and add constraint parameters with their unique names
      constraint_set.SetUniqueParameterNames();
      for (const auto& label : constraint_set.GetParameterNames())
      {
        if (params_map.count(label) > 0)
          throw MicmException(
              MicmSeverity::Error,
              MICM_ERROR_CATEGORY_SOLVER,
              MICM_SOLVER_ERROR_CODE_DUPLICATE_PARAMETER,
              "Duplicate parameter name: " + label);
        params_map.emplace(label, params_map.size());
      }

      algebraic_variable_ids = constraint_set.AlgebraicVariableIds();
      rates.SetAlgebraicVariableIds(algebraic_variable_ids);
      for (const auto variable_id : algebraic_variable_ids)
        mass_matrix_diagonal[variable_id] = 0.0;

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

    // Resolve external model constraints (runtime activation)
    if (!external_constraint_models_.empty())
    {
      auto ext_constraint_models_copy = external_constraint_models_;
      constraint_set.SetExternalConstraintModels(std::move(ext_constraint_models_copy));
      constraint_set.ResolveExternalConstraints(species_map);

      // Add external constraint parameter names to the params map
      for (const auto& label : constraint_set.ExternalConstraintParameterNames())
      {
        if (params_map.count(label) > 0)
          throw MicmException(
              MicmSeverity::Error,
              MICM_ERROR_CATEGORY_SOLVER,
              MICM_SOLVER_ERROR_CODE_DUPLICATE_PARAMETER,
              "Duplicate parameter name: " + label);
        params_map.emplace(label, params_map.size());
      }

      // Add initialize constraint parameter names to the params map
      for (const auto& label : constraint_set.ExternalInitializeConstraintParameterNames())
      {
        if (params_map.count(label) > 0)
          throw MicmException(
              MicmSeverity::Error,
              MICM_ERROR_CATEGORY_SOLVER,
              MICM_SOLVER_ERROR_CODE_DUPLICATE_PARAMETER,
              "Duplicate parameter name: " + label);
        params_map.emplace(label, params_map.size());
      }

      auto ext_algebraic_ids = constraint_set.AlgebraicVariableIds();
      // Find newly added algebraic IDs from external models
      for (const auto& id : ext_algebraic_ids)
      {
        if (algebraic_variable_ids.insert(id).second)
        {
          // Filter kinetic sparsity entries from this new algebraic row
          for (auto it = nonzero_elements.begin(); it != nonzero_elements.end();)
          {
            if (it->first == id)
              it = nonzero_elements.erase(it);
            else
              ++it;
          }
          mass_matrix_diagonal[id] = 0.0;
        }
      }
      rates.SetAlgebraicVariableIds(algebraic_variable_ids);

      // Merge external constraint Jacobian sparsity
      auto ext_jac_elements = constraint_set.ExternalNonZeroJacobianElements(species_map);
      nonzero_elements.insert(ext_jac_elements.begin(), ext_jac_elements.end());
    }

    // Re-add external model process Jacobian elements for algebraic rows.
    // Built-in ProcessSet is protected by is_algebraic_variable_ guards that skip
    // algebraic rows, but external models' JacobianFunction closures pre-compute
    // VectorIndex at setup time and need these elements to exist in the sparse matrix.
    if (!algebraic_variable_ids.empty())
    {
      for (const auto& model : external_models_)
      {
        auto ext_process_elements = model.non_zero_jacobian_elements_func_(species_map);
        for (const auto& elem : ext_process_elements)
        {
          if (algebraic_variable_ids.count(elem.first) > 0)
            nonzero_elements.insert(elem);
        }
      }
    }

    auto jacobian = BuildJacobian<SparseMatrixPolicy>(nonzero_elements, 1, number_of_species, true);

    LinearSolverPolicy linear_solver(jacobian, 0);
    if constexpr (LuDecompositionInPlaceConcept<LuDecompositionPolicy, SparseMatrixPolicy>)
    {
      auto lu = LuDecompositionPolicy::template GetLUMatrix<SparseMatrixPolicy>(jacobian, 0, true);
      jacobian = std::move(lu);
    }

    std::vector<std::string> variable_names{ number_of_species };
    for (auto& species_pair : species_map)
      variable_names[species_pair.second] = species_pair.first;

    // Build the params map after the constraint set is created,
    // since it adds its parameters to the map.
    std::vector<std::string> labels{ params_map.size() };
    for (auto& param_pair : params_map)
      labels[param_pair.second] = param_pair.first;

    rates.SetJacobianFlatIds(jacobian);
    rates.SetExternalModelFunctions(params_map, species_map, jacobian);

    // Compile external model update functions now that params_map is finalized
    for (const auto& model : external_models_)
    {
      update_state_param_funcs.push_back(model.update_state_parameters_function_(params_map));
    }

    std::vector<std::function<void(const DenseMatrixPolicy&, DenseMatrixPolicy&)>> init_constraint_param_funcs;

    if (constraint_set.Size() > 0)
    {
      constraint_set.SetJacobianFlatIds(jacobian);

      // Set forcing, jacobian, updating state param functions
      // The species map and parameter map are used to set indices in the state variables
      // and custom parameters.
      constraint_set.SetConstraintFunctions(species_map, params_map, jacobian);
      constraint_set.SetExternalModelConstraintFunctions(params_map, species_map, jacobian);

      // Add functions that update state parameters when temperature changes
      auto constraint_param_funcs = constraint_set.GetUpdateStateParamFunctions();
      update_state_param_funcs.insert(
          update_state_param_funcs.end(), constraint_param_funcs.begin(), constraint_param_funcs.end());

      // Collect constraint parameter initialization functions
      auto ext_init_funcs = constraint_set.GetExternalInitializeConstraintParamFunctions();
      init_constraint_param_funcs.insert(
          init_constraint_param_funcs.end(), ext_init_funcs.begin(), ext_init_funcs.end());
          
      // Add external constraint parameter update functions to the pipeline
      auto ext_constraint_param_funcs = constraint_set.GetExternalUpdateStateParamFunctions();
      update_state_param_funcs.insert(
          update_state_param_funcs.end(), ext_constraint_param_funcs.begin(), ext_constraint_param_funcs.end());
    }

    StateParameters state_parameters = { .number_of_species_ = number_of_species,
                                         .number_of_constraints_ = constraint_set.Size(),
                                         .number_of_rate_constants_ = this->reactions_.size(),
                                         .variable_names_ = variable_names,
                                         .custom_rate_parameter_labels_ = labels,
                                         .nonzero_jacobian_elements_ = nonzero_elements,
                                         .mass_matrix_diagonal_ = mass_matrix_diagonal };

    this->SetAbsoluteTolerances(state_parameters.absolute_tolerance_, species_map);

    // make a copy of the options so that the builder can be used repeatedly
    // this matters because the absolute tolerances must be set to match the system size, and that may change
    auto options = this->options_;

    return Solver<SolverPolicy, StatePolicy>(
        SolverPolicy(std::move(linear_solver), std::move(rates), std::move(constraint_set)),
        state_parameters,
        options,
        reactions_,
        system_,
        update_state_param_funcs,
        init_constraint_param_funcs);
  }

}  // namespace micm
