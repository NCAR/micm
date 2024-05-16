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
  inline SolverBuilder& SolverBuilder::SetSystem(micm::System system)
  {
    system_ = system;
    return *this;
  }

  inline SolverBuilder& SolverBuilder::SetReactions(std::vector<micm::Process> reactions)
  {
    reactions_ = reactions;
    return *this;
  }

  inline SolverBuilder& SolverBuilder::SetNumberOfGridCells(int number_of_grid_cells)
  {
    number_of_grid_cells_ = number_of_grid_cells;
    return *this;
  }

  inline SolverBuilder& SolverBuilder::SetSolverParameters(const RosenbrockSolverParameters& options)
  {
    if (!std::holds_alternative<std::monostate>(options_))
    {
      throw std::runtime_error("Solver type already set");
    }
    options_.emplace<RosenbrockSolverParameters>(options);
    return *this;
  }

  inline SolverBuilder& SolverBuilder::SetSolverParameters(const BackwardEulerSolverParameters& options)
  {
    if (!std::holds_alternative<std::monostate>(options_))
    {
      throw std::runtime_error("Solver type already set");
    }
    options_.emplace<BackwardEulerSolverParameters>(options);
    return *this;
  }

  inline Solver SolverBuilder::Build()
  {
    if (std::holds_alternative<micm::RosenbrockSolverParameters>(options_))
    {
      throw std::runtime_error("Not implemented yet");
    }
    if (std::holds_alternative<micm::BackwardEulerSolverParameters>(options_))
    {
      return BuildBackwardEulerSolver();
    }

    throw std::runtime_error("No solver type set");
  }

  template<class ProcessSetPolicy>
  inline void SolverBuilder::UnusedSpeciesCheck()
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

  template<class MatrixPolicy, class ProcessSetPolicy>
  inline std::map<std::string, std::size_t> SolverBuilder::GetSpeciesMap() const
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

      using Matrix = typename MatrixPolicy::IntMatrix;
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

  inline void SolverBuilder::SetAbsoluteTolerances(
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

  inline std::vector<std::string> SolverBuilder::GetCustomParameterLabels() const
  {
    std::vector<std::string> param_labels{};
    for (const auto& reaction : reactions_)
      if (reaction.rate_constant_)
        for (auto& label : reaction.rate_constant_->CustomParameters())
          param_labels.push_back(label);
    return param_labels;
  }

  inline std::vector<std::size_t> SolverBuilder::GetJacobianDiagonalElements(auto jacobian) const
  {
    std::vector<std::size_t> jacobian_diagonal_elements;

    jacobian_diagonal_elements.reserve(jacobian.NumRows());

    for (std::size_t i = 0; i < jacobian.NumRows(); ++i)
    {
      jacobian_diagonal_elements.push_back(jacobian.VectorIndex(0, i, i));
    }

    return jacobian_diagonal_elements;
  }

  template<class MatrixPolicy, class SparseMatrixPolicy>
  inline Solver CpuSolverBuilder<MatrixPolicy, SparseMatrixPolicy>::BuildBackwardEulerSolver()
  {
    using ProcessSetPolicy = micm::ProcessSet;

    auto parameters = std::get<BackwardEulerSolverParameters>(options_);
    auto species_map = GetSpeciesMap<MatrixPolicy, ProcessSetPolicy>();
    auto labels = GetCustomParameterLabels();
    std::size_t number_of_species = system_.StateSize();

    UnusedSpeciesCheck<ProcessSetPolicy>();
    SetAbsoluteTolerances(parameters.absolute_tolerance_, species_map);

    ProcessSetPolicy process_set(reactions_, species_map);
    auto jacobian =
        BuildJacobian<SparseMatrixPolicy>(process_set.NonZeroJacobianElements(), number_of_grid_cells_, number_of_species);
    auto diagonal_elements = GetJacobianDiagonalElements(jacobian);

    process_set.SetJacobianFlatIds(jacobian);
    micm::LinearSolver<typename MatrixPolicy::value_type, SparseMatrixPolicy> linear_solver(jacobian, 1e-30);

    return Solver(
        new SolverImpl<decltype(linear_solver), decltype(process_set)>(),
        parameters,
        number_of_grid_cells_,
        number_of_species,
        reactions_.size());
  }

}  // namespace micm