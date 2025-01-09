// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

enum class MicmStateErrc
{
  UnknownSpecies = 1,                                                   // Unknown species
  UnknownRateConstantParameter = 2,                                     // Unknown rate constant parameter
  IncorrectNumberOfConcentrationValuesForMultiGridcellState = 3,        // Incorrect number of concentration values
  IncorrectNumberOfCustomRateParameterValues = 4,                       // Incorrect number of custom rate parameter values
  IncorrectNumberOfCustomRateParameterValuesForMultiGridcellState = 5,  // Incorrect number of grid cells
};

namespace std
{
  template<>
  struct is_error_condition_enum<MicmStateErrc> : true_type
  {
  };
}  // namespace std

namespace
{

  class MicmStateErrorCategory : public std::error_category
  {
   public:
    const char* name() const noexcept override
    {
      return "MICM State";
    }
    std::string message(int ev) const override
    {
      switch (static_cast<MicmStateErrc>(ev))
      {
        case MicmStateErrc::UnknownSpecies: return "Unknown species";
        case MicmStateErrc::UnknownRateConstantParameter: return "Unknown rate constant parameter";
        case MicmStateErrc::IncorrectNumberOfConcentrationValuesForMultiGridcellState:
          return "Incorrect number of concentration values for multi-gridcell State";
        case MicmStateErrc::IncorrectNumberOfCustomRateParameterValues:
          return "Incorrect number of custom rate parameter values per grid cell";
        case MicmStateErrc::IncorrectNumberOfCustomRateParameterValuesForMultiGridcellState:
          return "Incorrect number of custom rate parameter values for multi-gridcell State";
        default: return "Unknown error";
      }
    }
  };

  const MicmStateErrorCategory micmStateErrorCategory{};

}  // namespace

inline std::error_code make_error_code(MicmStateErrc e)
{
  return { static_cast<int>(e), micmStateErrorCategory };
}

namespace micm
{

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
  inline State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>::State()
      : conditions_(),
        variables_(),
        custom_rate_parameters_(),
        rate_constants_(),
        jacobian_()
  {
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
  inline State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>::State(
      const StateParameters& parameters)
      : conditions_(parameters.number_of_grid_cells_),
        variables_(parameters.number_of_grid_cells_, parameters.variable_names_.size(), 0.0),
        custom_rate_parameters_(parameters.number_of_grid_cells_, parameters.custom_rate_parameter_labels_.size(), 0.0),
        rate_constants_(parameters.number_of_grid_cells_, parameters.number_of_rate_constants_, 0.0),
        variable_map_(),
        custom_rate_parameter_map_(),
        variable_names_(parameters.variable_names_),
        jacobian_(),
        lower_matrix_(),
        upper_matrix_(),
        state_size_(parameters.variable_names_.size()),
        number_of_grid_cells_(parameters.number_of_grid_cells_)
  {
    std::size_t index = 0;
    for (auto& name : variable_names_)
      variable_map_[name] = index++;
    index = 0;
    for (auto& label : parameters.custom_rate_parameter_labels_)
      custom_rate_parameter_map_[label] = index++;

    jacobian_ = BuildJacobian<SparseMatrixPolicy>(
        parameters.nonzero_jacobian_elements_, parameters.number_of_grid_cells_, state_size_);
    
    if constexpr (LuDecompositionInPlaceConcept<LuDecompositionPolicy, SparseMatrixPolicy>)
    {
      auto lu = LuDecompositionPolicy::template GetLUMatrix<SparseMatrixPolicy>(jacobian_, 0);
      jacobian_ = std::move(lu);
    }
    else
    {
      auto lu = LuDecompositionPolicy::template GetLUMatrices<SparseMatrixPolicy, LMatrixPolicy, UMatrixPolicy>(jacobian_, 0);
      auto lower_matrix = std::move(lu.first);
      auto upper_matrix = std::move(lu.second);
      lower_matrix_ = lower_matrix;
      upper_matrix_ = upper_matrix;
    }
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
  inline void
  State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>::SetConcentrations(
      const std::unordered_map<std::string, std::vector<double>>& species_to_concentration)
  {
    const std::size_t num_grid_cells = conditions_.size();
    for (const auto& pair : species_to_concentration)
      SetConcentration({ pair.first }, pair.second);
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
  inline void
  State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>::SetConcentration(
      const Species& species,
      double concentration)
  {
    auto var = variable_map_.find(species.name_);
    if (var == variable_map_.end())
      throw std::system_error(make_error_code(MicmStateErrc::UnknownSpecies), species.name_);
    if (variables_.NumRows() != 1)
      throw std::system_error(make_error_code(MicmStateErrc::IncorrectNumberOfConcentrationValuesForMultiGridcellState));
    variables_[0][variable_map_[species.name_]] = concentration;
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
  inline void
  State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>::SetConcentration(
      const Species& species,
      const std::vector<double>& concentration)
  {
    auto var = variable_map_.find(species.name_);
    if (var == variable_map_.end())
      throw std::system_error(make_error_code(MicmStateErrc::UnknownSpecies), species.name_);
    if (variables_.NumRows() != concentration.size())
      throw std::system_error(make_error_code(MicmStateErrc::IncorrectNumberOfConcentrationValuesForMultiGridcellState));
    std::size_t i_species = variable_map_[species.name_];
    for (std::size_t i = 0; i < variables_.NumRows(); ++i)
      variables_[i][i_species] = concentration[i];
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
  inline void State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>::
      UnsafelySetCustomRateParameters(const std::vector<std::vector<double>>& parameters)
  {
    if (parameters.size() != variables_.NumRows())
      throw std::system_error(
          make_error_code(MicmStateErrc::IncorrectNumberOfCustomRateParameterValuesForMultiGridcellState));

    if (parameters[0].size() != custom_rate_parameters_.NumColumns())
      throw std::system_error(make_error_code(MicmStateErrc::IncorrectNumberOfCustomRateParameterValues));

    for (size_t i = 0; i < number_of_grid_cells_; ++i)
    {
      custom_rate_parameters_[i] = parameters[i];
    }
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
  inline void
  State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>::SetCustomRateParameters(
      const std::unordered_map<std::string, std::vector<double>>& parameters)
  {
    for (auto& pair : parameters)
      SetCustomRateParameter(pair.first, pair.second);
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
  inline void
  State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>::SetCustomRateParameter(
      const std::string& label,
      double value)
  {
    auto param = custom_rate_parameter_map_.find(label);
    if (param == custom_rate_parameter_map_.end())
      throw std::system_error(make_error_code(MicmStateErrc::UnknownRateConstantParameter), label);
    if (custom_rate_parameters_.NumRows() != 1)
      throw std::system_error(
          make_error_code(MicmStateErrc::IncorrectNumberOfCustomRateParameterValuesForMultiGridcellState));
    custom_rate_parameters_[0][param->second] = value;
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
  inline void
  State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>::SetCustomRateParameter(
      const std::string& label,
      const std::vector<double>& values)
  {
    auto param = custom_rate_parameter_map_.find(label);
    if (param == custom_rate_parameter_map_.end())
      throw std::system_error(make_error_code(MicmStateErrc::UnknownRateConstantParameter), label);
    if (custom_rate_parameters_.NumRows() != values.size())
      throw std::system_error(
          make_error_code(MicmStateErrc::IncorrectNumberOfCustomRateParameterValuesForMultiGridcellState));
    for (std::size_t i = 0; i < custom_rate_parameters_.NumRows(); ++i)
      custom_rate_parameters_[i][param->second] = values[i];
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
  inline void
  State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>::PrintHeader()
  {
    auto largest_str_iter = std::max_element(
        variable_names_.begin(), variable_names_.end(), [](const auto& a, const auto& b) { return a.size() < b.size(); });
    int largest_str_size = largest_str_iter->size();
    int width = (largest_str_size < 10) ? 11 : largest_str_size + 2;

    std::cout << std::setw(6) << "time";
    if (variables_.NumRows() > 1)
    {
      std::cout << "," << std::setw(6) << "grid";
    }

    for (const auto& [species, index] : variable_map_)
    {
      std::cout << "," << std::setw(width) << species;
    }
    std::cout << std::endl;
  }

  template<
      class DenseMatrixPolicy,
      class SparseMatrixPolicy,
      class LuDecompositionPolicy,
      class LMatrixPolicy,
      class UMatrixPolicy>
  inline void State<DenseMatrixPolicy, SparseMatrixPolicy, LuDecompositionPolicy, LMatrixPolicy, UMatrixPolicy>::PrintState(
      double time)
  {
    std::ios oldState(nullptr);
    oldState.copyfmt(std::cout);

    auto largest_str_iter = std::max_element(
        variable_names_.begin(), variable_names_.end(), [](const auto& a, const auto& b) { return a.size() < b.size(); });
    int largest_str_size = largest_str_iter->size();
    int width = (largest_str_size < 10) ? 11 : largest_str_size + 2;

    for (size_t i = 0; i < variables_.NumRows(); ++i)
    {
      std::cout << std::setw(6) << time << ",";

      if (variables_.NumRows() > 1)
      {
        std::cout << std::setw(6) << i << ",";
      }

      bool first = true;
      for (const auto& [species, index] : variable_map_)
      {
        if (!first)
        {
          std::cout << ",";
        }
        std::cout << std::scientific << std::setw(width) << std::setprecision(2) << variables_[i][index];
        first = false;
      }
      std::cout << std::endl;
      std::cout.copyfmt(oldState);
    }
  }
}  // namespace micm
