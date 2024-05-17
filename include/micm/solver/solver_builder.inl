/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */

enum class MicmSolverBuilderErrc
{
  UnusedSpecies = 1,  // Unused species present in the chemical system
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
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  inline SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>&
  SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::SetSystem(System system)
  {
    system_ = system;
    return *this;
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  inline SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>&
  SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::SetReactions(std::vector<Process> reactions)
  {
    reactions_ = reactions;
    return *this;
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  inline SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>&
  SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::SetNumberOfGridCells(int number_of_grid_cells)
  {
    number_of_grid_cells_ = number_of_grid_cells;
    return *this;
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  inline SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>&
  SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::SetSolverParameters(const RosenbrockSolverParameters& options)
  {
    if (!std::holds_alternative<std::monostate>(options_))
    {
      throw std::runtime_error("Solver type already set");
    }
    options_.emplace<RosenbrockSolverParameters>(options);
    return *this;
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  inline SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>&
  SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::SetSolverParameters(const BackwardEulerSolverParameters& options)
  {
    if (!std::holds_alternative<std::monostate>(options_))
    {
      throw std::runtime_error("Solver type already set");
    }
    options_.emplace<BackwardEulerSolverParameters>(options);
    return *this;
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  inline auto SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::Build()
  {
    if (std::holds_alternative<RosenbrockSolverParameters>(options_))
    {
      throw std::runtime_error("Not implemented yet");
    }
    if (std::holds_alternative<BackwardEulerSolverParameters>(options_))
    {
      return BuildBackwardEulerSolver();
    }

    throw std::runtime_error("No solver type set");
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  template<class ProcessSetPolicy>
  inline void SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::UnusedSpeciesCheck()
  {
    if (ignore_unused_species_)
    {
      return;
    }

    std::size_t index = 0;
    auto used_species = ProcessSetPolicy::SpeciesUsed(reactions_);
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

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  template<class ProcessSetPolicy>
  inline std::map<std::string, std::size_t> SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::GetSpeciesMap() const
  {
    std::map<std::string, std::size_t> species_map;
    std::function<std::string(const std::vector<std::string>& variables, const std::size_t i)> state_reordering;
    std::size_t index = 0;
    for (auto& name : system_.UniqueNames())
      species_map[name] = index++;

    if (reorder_state_)
    {
      // get unsorted Jacobian non-zero elements
      auto unsorted_process_set = ProcessSetPolicy(reactions_, species_map);
      auto unsorted_jac_elements = unsorted_process_set.NonZeroJacobianElements();

      using Matrix = typename DenseMatrixPolicy::IntMatrix;
      Matrix unsorted_jac_non_zeros(system_.StateSize(), system_.StateSize(), 0);
      for (auto& elem : unsorted_jac_elements)
        unsorted_jac_non_zeros[elem.first][elem.second] = 1;
      auto reorder_map = DiagonalMarkowitzReorder<Matrix>(unsorted_jac_non_zeros);

      state_reordering = [=](const std::vector<std::string>& variables, const std::size_t i)
      { return variables[reorder_map[i]]; };

      std::size_t index = 0;
      for (auto& name : system_.UniqueNames(state_reordering))
        species_map[name] = index++;
    }

    return species_map;
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  inline void SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::SetAbsoluteTolerances(
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
          tolerances[species_map.at(species.name_)] = species.GetProperty<double>("absolute tolerance");
        }
      }
      for (auto& phase : system_.phases_)
      {
        for (auto& species : phase.second.species_)
        {
          if (species.HasProperty("absolute tolerance"))
          {
            tolerances[species_map.at(species.name_)] = species.GetProperty<double>("absolute tolerance");
          }
        }
      }
    }
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  inline std::vector<std::string> SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::GetCustomParameterLabels() const
  {
    std::vector<std::string> param_labels{};
    for (const auto& reaction : reactions_)
      if (reaction.rate_constant_)
        for (auto& label : reaction.rate_constant_->CustomParameters())
          param_labels.push_back(label);
    return param_labels;
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  inline std::vector<std::size_t> SolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::GetJacobianDiagonalElements(
      auto jacobian) const
  {
    std::vector<std::size_t> jacobian_diagonal_elements;

    jacobian_diagonal_elements.reserve(jacobian.NumRows());

    for (std::size_t i = 0; i < jacobian.NumRows(); ++i)
    {
      jacobian_diagonal_elements.push_back(jacobian.VectorIndex(0, i, i));
    }

    return jacobian_diagonal_elements;
  }

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  inline Solver<
      BackwardEuler<LinearSolver<typename DenseMatrixPolicy::value_type, SparseMatrixPolicy>, ProcessSet>,
      State<DenseMatrixPolicy, SparseMatrixPolicy>>
  CpuSolverBuilder<DenseMatrixPolicy, SparseMatrixPolicy>::BuildBackwardEulerSolver()
  {
    using ProcessSetPolicy = ProcessSet;
    using LinearSolverPolicy = LinearSolver<typename DenseMatrixPolicy::value_type, SparseMatrixPolicy>;

    auto parameters = std::get<BackwardEulerSolverParameters>(this->options_);
    auto species_map = this->template GetSpeciesMap<ProcessSetPolicy>();
    auto labels = this->GetCustomParameterLabels();
    std::size_t number_of_species = this->system_.StateSize();

    this->template UnusedSpeciesCheck<ProcessSetPolicy>();
    this->SetAbsoluteTolerances(parameters.absolute_tolerance_, species_map);

    ProcessSetPolicy process_set(this->reactions_, species_map);
    auto diagonal_elements = process_set.NonZeroJacobianElements();
    auto jacobian = BuildJacobian<SparseMatrixPolicy>(diagonal_elements, this->number_of_grid_cells_, number_of_species);
    auto jacobian_diagonal_elements = this->GetJacobianDiagonalElements(jacobian);

    process_set.SetJacobianFlatIds(jacobian);
    LinearSolverPolicy linear_solver(jacobian, 1e-30);

    std::vector<std::string> variable_names{ number_of_species };
    for (auto& species_pair : species_map)
      variable_names[species_pair.second] = species_pair.first;

    StateParameters state_parameters_ = { .number_of_grid_cells_ = this->number_of_grid_cells_,
                                          .number_of_species_ = number_of_species,
                                          .number_of_rate_constants_ = this->reactions_.size(),
                                          .variable_names_ = variable_names,
                                          .custom_rate_parameter_labels_ = labels,
                                          .jacobian_diagonal_elements_ = jacobian_diagonal_elements,
                                          .nonzero_jacobian_elements_ = diagonal_elements };

    return Solver<BackwardEuler<LinearSolverPolicy, ProcessSetPolicy>, State<DenseMatrixPolicy, SparseMatrixPolicy>>(
        BackwardEuler<LinearSolverPolicy, ProcessSetPolicy>(
            parameters, linear_solver, process_set, jacobian_diagonal_elements, this->reactions_),
        state_parameters_,
        this->number_of_grid_cells_,
        number_of_species,
        this->reactions_.size());
  }

}  // namespace micm