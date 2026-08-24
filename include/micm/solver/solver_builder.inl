// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/solver/external_model_dispatcher.hpp>
#include <micm/util/types.hpp>

namespace micm
{
  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy,
      class... ExternalModels>
  inline void SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy,
      ExternalModels...>::UnusedSpeciesCheck(const RatesPolicy& rates) const
  {
    if (ignore_unused_species_)
    {
      return;
    }

    auto used_species = rates.SpeciesUsed(reactions_);
    for (const auto& constraint : constraints_)
    {
      for (const auto& dep : constraint.SpeciesDependencies())
      {
        used_species.insert(dep);
      }
      used_species.insert(constraint.AlgebraicSpecies());
    }
    for (const auto& ps : external_process_sets_)
    {
      auto ext = ps.species_used_func_();
      used_species.insert(ext.begin(), ext.end());
    }
    for (const auto& constraint : external_constraints_)
    {
      auto deps = constraint.species_dependencies_func_();
      used_species.insert(deps.begin(), deps.end());
      auto alg_names = constraint.algebraic_variable_names_func_();
      used_species.insert(alg_names.begin(), alg_names.end());
    }

    auto available_species = this->MergedUniqueNames();
    std::sort(available_species.begin(), available_species.end());
    std::set<std::string> unused_species;
    std::set_difference(
        available_species.begin(),
        available_species.end(),
        used_species.begin(),
        used_species.end(),
        std::inserter(unused_species, unused_species.begin()));
    if (!unused_species.empty())
    {
      std::string err_msg = "Unused species in chemical system:";
      for (const auto& species : unused_species)
      {
        err_msg += " '" + species + "'";
      }
      err_msg += ".";
      throw MicmException(MICM_ERROR_CATEGORY_SOLVER, MICM_SOLVER_ERROR_CODE_UNUSED_SPECIES, err_msg);
    }
  }

  template<
      class SolverParametersPolicy,
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class RatesPolicy,
      class LuDecompositionPolicy,
      class LinearSolverPolicy,
      class StatePolicy,
      class... ExternalModels>
  inline std::unordered_map<std::string, Index> SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy,
      ExternalModels...>::GetSpeciesMap() const
  {
    std::unordered_map<std::string, Index> species_map;
    Index index = 0;

    auto all_names = this->MergedUniqueNames();
    for (auto& name : all_names)
    {
      species_map[name] = index++;
    }

    if (reorder_state_)
    {
      auto unsorted_rates = RatesPolicy(reactions_, species_map);
      auto unsorted_jac_elements = unsorted_rates.NonZeroJacobianElements();
      for (const auto& ps : external_process_sets_)
      {
        auto ext = ps.non_zero_jacobian_elements_func_(species_map);
        unsorted_jac_elements.insert(ext.begin(), ext.end());
      }

      using Matrix = typename DenseMatrixPolicy::IntMatrix;
      const auto n = this->MergedStateSize();
      Matrix unsorted_jac_non_zeros(n, n, 0);
      for (auto& elem : unsorted_jac_elements)
      {
        unsorted_jac_non_zeros[elem.first][elem.second] = 1;
      }
      auto reorder_map = DiagonalMarkowitzReorder<Matrix>(unsorted_jac_non_zeros);

      index = 0;
      for (Index i = 0; i < all_names.size(); ++i)
      {
        species_map[all_names[reorder_map[i]]] = index++;
      }
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
      class StatePolicy,
      class... ExternalModels>
  inline std::unordered_map<std::string, Index> SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy,
      ExternalModels...>::GetCustomParameterMap() const
  {
    std::unordered_map<std::string, Index> params{};
    std::vector<std::string> duplicates;

    auto add_param = [&params, &duplicates](const std::string& label, const std::string& source)
    {
      auto [it, added] = params.emplace(label, params.size());
      if (!added)
      {
        duplicates.push_back(label + " (from " + source + ")");
      }
    };

    for (const auto& reaction : reactions_)
    {
      const auto& process = reaction.process_;
      if (const auto* ud = std::get_if<UserDefinedRateConstantParameters>(&process.rate_constant_))
      {
        add_param(ud->label_, "reaction");
      }
      else if (const auto* surf = std::get_if<SurfaceRateConstantParameters>(&process.rate_constant_))
      {
        add_param(surf->label_ + ".effective radius [m]", "reaction");
        add_param(surf->label_ + ".particle number concentration [# m-3]", "reaction");
      }
    }

    for (const auto& sys : external_systems_)
    {
      auto param_names = sys.parameter_names_func_();
      for (const auto& label : param_names)
      {
        add_param(label, "external_model");
      }
    }

    if (!duplicates.empty())
    {
      std::ostringstream oss;
      oss << "Duplicate parameter labels detected:\n";
      for (const auto& d : duplicates)
      {
        oss << "  - " << d << "\n";
      }

      throw MicmException(MICM_ERROR_CATEGORY_SOLVER, MICM_SOLVER_ERROR_CODE_DUPLICATE_PARAMETER, oss.str());
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
      class StatePolicy,
      class... ExternalModels>
  inline void SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy,
      ExternalModels...>::
      SetAbsoluteTolerances(std::vector<Real>& tolerances, const std::unordered_map<std::string, Index>& species_map) const
  {
    tolerances = std::vector<Real>(species_map.size(), 1e-3);
    for (const auto& phase_species : system_.gas_phase_.phase_species_)
    {
      const auto& species = phase_species.species_;
      if (species.HasProperty("absolute tolerance"))
      {
        tolerances[species_map.at(species.name_)] = species.template GetProperty<Real>("absolute tolerance");
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
      class StatePolicy,
      class... ExternalModels>
  inline auto SolverBuilder<
      SolverParametersPolicy,
      DenseMatrixPolicy,
      SparseMatrixPolicy,
      RatesPolicy,
      LuDecompositionPolicy,
      LinearSolverPolicy,
      StatePolicy,
      ExternalModels...>::Build()
  {
    if (!valid_system_)
    {
      throw MicmException(
          MICM_ERROR_CATEGORY_SOLVER, MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SYSTEM, "Missing chemical system.");
    }

    Index number_of_species = this->MergedStateSize();
    if (number_of_species == 0)
    {
      throw MicmException(
          MICM_ERROR_CATEGORY_SOLVER,
          MICM_SOLVER_ERROR_CODE_MISSING_CHEMICAL_SPECIES,
          "Provided chemical system contains no species.");
    }

    // Constraints are not supported with CUDA matrix policies
    constexpr bool is_cuda_policy = requires(DenseMatrixPolicy m) {
      m.CopyToDevice();
      m.CopyToHost();
      m.AsDeviceParam();
    };
    if constexpr (is_cuda_policy)
    {
      if (!constraints_.empty() || !external_constraints_.empty())
      {
        throw MicmException(
            MICM_ERROR_CATEGORY_SOLVER,
            MICM_SOLVER_ERROR_CODE_CUDA_CONSTRAINTS_UNSUPPORTED,
            "Constraints are not supported with CUDA matrix policies.");
      }
    }

    using InnerConstraintSet = ConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>;
    using RatesBundleType = RatesBundle<RatesPolicy, ExternalModels...>;
    using ConstraintBundleType = ConstraintBundle<InnerConstraintSet, ExternalModels...>;
    using SolverPolicy =
        typename SolverParametersPolicy::template SolverType<RatesBundleType, LinearSolverPolicy, ConstraintBundleType>;

    // Sort reactions by rate constant type so ReactionRateConstantStore and ProcessSet
    // share a consistent ordering and each type occupies a contiguous block.
    std::stable_sort(
        reactions_.begin(),
        reactions_.end(),
        [](const Process& a, const Process& b)
        {
          auto type_order = [](const Process& p) -> int
          {
            return std::visit(
                [](const auto& v) -> int
                {
                  using T = std::decay_t<decltype(v)>;
                  if constexpr (std::is_same_v<T, ArrheniusRateConstantParameters>)
                  {
                    return 0;
                  }
                  else if constexpr (std::is_same_v<T, TroeRateConstantParameters>)
                  {
                    return 1;
                  }
                  else if constexpr (std::is_same_v<T, TernaryChemicalActivationRateConstantParameters>)
                  {
                    return 2;
                  }
                  else if constexpr (std::is_same_v<T, BranchedRateConstantParameters>)
                  {
                    return 3;
                  }
                  else if constexpr (std::is_same_v<T, TunnelingRateConstantParameters>)
                  {
                    return 4;
                  }
                  else if constexpr (std::is_same_v<T, TaylorSeriesRateConstantParameters>)
                  {
                    return 5;
                  }
                  else if constexpr (std::is_same_v<T, ReversibleRateConstantParameters>)
                  {
                    return 6;
                  }
                  else if constexpr (std::is_same_v<T, UserDefinedRateConstantParameters>)
                  {
                    return 7;
                  }
                  else if constexpr (std::is_same_v<T, SurfaceRateConstantParameters>)
                  {
                    return 8;
                  }
                  else
                  {
                    return 9;  // LambdaRateConstantParameters
                  }
                },
                p.process_.rate_constant_);
          };
          return type_order(a) < type_order(b);
        });

    auto species_map = this->GetSpeciesMap();
    RatesPolicy rates(reactions_, species_map);
    this->UnusedSpeciesCheck(rates);

    // Build the ODE Jacobian sparsity, including contributions from external process models.
    auto nonzero_elements = rates.NonZeroJacobianElements();
    for (const auto& ps : external_process_sets_)
    {
      auto ext = ps.non_zero_jacobian_elements_func_(species_map);
      nonzero_elements.insert(ext.begin(), ext.end());
    }

    auto params_map = this->GetCustomParameterMap();

    InnerConstraintSet constraint_set;

    // Build mass-matrix diagonal: species rows default to ODE (1), rows replaced by constraints are algebraic (0).
    std::vector<Real> mass_matrix_diagonal(number_of_species, 1.0);
    std::set<Index> algebraic_variable_ids;

    if (!constraints_.empty())
    {
      constraint_set = InnerConstraintSet(std::move(constraints_), species_map);

      constraint_set.SetUniqueParameterNames();
      for (const auto& label : constraint_set.GetParameterNames())
      {
        if (params_map.count(label) > 0)
        {
          throw MicmException(
              MICM_ERROR_CATEGORY_SOLVER, MICM_SOLVER_ERROR_CODE_DUPLICATE_PARAMETER, "Duplicate parameter name: " + label);
        }
        params_map.emplace(label, params_map.size());
      }

      algebraic_variable_ids = constraint_set.AlgebraicVariableIds();
      for (const auto variable_id : algebraic_variable_ids)
      {
        mass_matrix_diagonal[variable_id] = 0.0;
      }

      for (auto it = nonzero_elements.begin(); it != nonzero_elements.end();)
      {
        if (algebraic_variable_ids.count(it->first) > 0)
        {
          it = nonzero_elements.erase(it);
        }
        else
        {
          ++it;
        }
      }

      auto constraint_jac_elements = constraint_set.NonZeroJacobianElements();
      nonzero_elements.insert(constraint_jac_elements.begin(), constraint_jac_elements.end());
    }

    // Resolve external constraint models: determine which are active and collect their contributions.
    std::array<bool, sizeof...(ExternalModels)> constraint_active_mask{};
    std::set<Index> external_algebraic_variable_ids;
    if (!external_constraints_.empty())
    {
      std::set<std::string> seen_param_names;
      std::set<std::string> seen_init_names;
      for (Index i = 0; i < external_constraints_.size(); ++i)
      {
        const auto& model = external_constraints_[i];
        auto alg_names = model.algebraic_variable_names_func_();
        if (alg_names.empty())
        {
          continue;
        }
        constraint_active_mask[i] = true;

        for (const auto& name : alg_names)
        {
          auto it = species_map.find(name);
          if (it == species_map.end())
          {
            throw MicmException(
                MICM_ERROR_CATEGORY_CONSTRAINT,
                MICM_CONSTRAINT_ERROR_CODE_UNKNOWN_SPECIES,
                "External model constraint targets unknown algebraic species '" + name + "'");
          }
          if (algebraic_variable_ids.count(it->second) > 0)
          {
            throw MicmException(
                MICM_ERROR_CATEGORY_CONSTRAINT,
                MICM_CONSTRAINT_ERROR_CODE_DUPLICATE_ALGEBRAIC_SPECIES,
                "Multiple constraints map to the same algebraic species row '" + name + "'");
          }
          external_algebraic_variable_ids.insert(it->second);
          algebraic_variable_ids.insert(it->second);
        }

        for (const auto& label : model.state_parameter_names_func_())
        {
          if (!seen_param_names.insert(label).second)
          {
            throw MicmException(
                MICM_ERROR_CATEGORY_CONSTRAINT,
                MICM_CONSTRAINT_ERROR_CODE_DUPLICATE_PARAMETER,
                "Duplicate external constraint parameter name across models: " + label);
          }
          if (params_map.count(label) > 0)
          {
            throw MicmException(
                MICM_ERROR_CATEGORY_SOLVER,
                MICM_SOLVER_ERROR_CODE_DUPLICATE_PARAMETER,
                "Duplicate parameter name: " + label);
          }
          params_map.emplace(label, params_map.size());
        }

        for (const auto& label : model.initialize_constraint_parameter_names_func_())
        {
          if (!seen_init_names.insert(label).second)
          {
            throw MicmException(
                MICM_ERROR_CATEGORY_CONSTRAINT,
                MICM_CONSTRAINT_ERROR_CODE_DUPLICATE_PARAMETER,
                "Duplicate external initialize constraint parameter name across models: " + label);
          }
          if (params_map.count(label) > 0)
          {
            throw MicmException(
                MICM_ERROR_CATEGORY_SOLVER,
                MICM_SOLVER_ERROR_CODE_DUPLICATE_PARAMETER,
                "Duplicate parameter name: " + label);
          }
          params_map.emplace(label, params_map.size());
        }

        auto ext_jac_elements = model.non_zero_jacobian_elements_func_(species_map);
        nonzero_elements.insert(ext_jac_elements.begin(), ext_jac_elements.end());
        for (const auto& name : alg_names)
        {
          auto row = species_map.at(name);
          nonzero_elements.insert(std::make_pair(row, row));
        }
      }

      for (const auto id : external_algebraic_variable_ids)
      {
        mass_matrix_diagonal[id] = 0.0;
        // Purge any kinetic sparsity remnants for these newly-algebraic rows.
        for (auto it = nonzero_elements.begin(); it != nonzero_elements.end();)
        {
          if (it->first == id && algebraic_variable_ids.count(it->second) == 0
              && it->second == id)  // keep diagonal
          {
            ++it;
          }
          else if (it->first == id)
          {
            // Row-column entry from external constraint: keep. Kinetic ODE contributions
            // that shared this row are dropped below because they were added from the
            // built-in ProcessSet (which now guards with is_algebraic_variable_).
            ++it;
          }
          else
          {
            ++it;
          }
        }
      }
    }

    // Push the merged algebraic-variable set into the inner rates + constraint set.
    rates.SetAlgebraicVariableIds(algebraic_variable_ids);
    if (!external_algebraic_variable_ids.empty())
    {
      constraint_set.AddExternalAlgebraicVariableIds(external_algebraic_variable_ids);
    }

    auto jacobian = BuildJacobian<SparseMatrixPolicy>(nonzero_elements, 1, number_of_species, true);

    LinearSolverPolicy linear_solver(jacobian, 0);
    if constexpr (LuDecompositionInPlaceConcept<LuDecompositionPolicy, SparseMatrixPolicy>)
    {
      auto lu = LuDecompositionPolicy::GetLUMatrix(jacobian, 0, true);
      jacobian = std::move(lu);
    }

    std::vector<std::string> variable_names(number_of_species);
    for (auto& species_pair : species_map)
    {
      variable_names[species_pair.second] = species_pair.first;
    }

    std::vector<std::string> labels(params_map.size());
    for (auto& param_pair : params_map)
    {
      labels[param_pair.second] = param_pair.first;
    }

    rates.SetJacobianFlatIds(jacobian);

    if (constraint_set.Size() > 0 || !external_algebraic_variable_ids.empty())
    {
      if (constraint_set.Size() > 0)
      {
        constraint_set.SetJacobianFlatIds(jacobian);
        constraint_set.SetConstraintFunctions(params_map);
      }
      constraint_set.FinalizeAlgebraicErrorFunction();
    }

    // Move concrete external models into shared ownership; both bundles refer to the same tuple.
    auto shared_models = std::make_shared<std::tuple<ExternalModels...>>(std::move(external_models_));

    // Give each participating model its build-time indices and Jacobian sparsity handles.
    std::apply(
        [&](auto&... m)
        {
          auto finalize = [&](auto& model)
          {
            using M = std::decay_t<decltype(model)>;
            if constexpr (HasProcesses<M>)
            {
              if constexpr (requires { model.FinalizeProcessSetup(params_map, species_map, jacobian); })
              {
                model.FinalizeProcessSetup(params_map, species_map, jacobian);
              }
            }
            if constexpr (HasConstraints<M>)
            {
              if constexpr (requires { model.FinalizeConstraintSetup(params_map, species_map, jacobian); })
              {
                model.FinalizeConstraintSetup(params_map, species_map, jacobian);
              }
            }
          };
          (finalize(m), ...);
        },
        *shared_models);

    RatesBundleType rates_bundle(std::move(rates), shared_models);
    ConstraintBundleType constraint_bundle(std::move(constraint_set), shared_models, constraint_active_mask);

    StateParameters state_parameters = { .number_of_species_ = number_of_species,
                                         .number_of_constraints_ = constraint_bundle.Size(),
                                         .number_of_rate_constants_ = this->reactions_.size(),
                                         .variable_names_ = variable_names,
                                         .custom_rate_parameter_labels_ = labels,
                                         .nonzero_jacobian_elements_ = nonzero_elements,
                                         .mass_matrix_diagonal_ = mass_matrix_diagonal };

    this->SetAbsoluteTolerances(state_parameters.absolute_tolerance_, species_map);

    auto options = this->options_;

    return Solver<SolverPolicy, StatePolicy>(
        SolverPolicy(std::move(linear_solver), std::move(rates_bundle), std::move(constraint_bundle)),
        state_parameters,
        options,
        reactions_,
        system_);
  }

}  // namespace micm
