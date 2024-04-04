#include <micm/configure/solver_config.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/profiler/instrumentation.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>
#include "initial_conditions.hpp"

#include <iomanip>
#include <cassert>
#include <iostream>
#include <random>
#include <tuple>
#include <chrono>
#include <sstream>
#include <filesystem>


namespace fs = std::filesystem;
using namespace micm;

void printMap(const std::unordered_map<std::string, std::vector<double>> &data_map)
{
  for (const auto &entry : data_map)
  {
    std::cout << "Key: " << entry.first << ", Values: ";
    for (const double value : entry.second)
    {
        std::cout << value << " ";
    }
    std::cout << std::endl;
  }
}

int OutputSpeciesState(const std::string file_prefix, const std::vector<std::tuple<std::string, double, double>> &species)
{
  fs::path file_path = file_prefix + "species_state.csv";
  std::ofstream csv_file(file_path);
  if (!csv_file.is_open())
  {
    std::cerr << "Failed to open the CSV file " << file_path << "." << std::endl;
    return 1;
  }
  csv_file << "Species Name,Initial Concentration,Initial Forcing" << std::endl;
  for (auto &species_data : species)
    csv_file << std::get<0>(species_data) << ","
            << std::setprecision(10)
            << std::get<1>(species_data) << ","
            << std::get<2>(species_data) << std::endl;
  csv_file.close();
  return 0;
}

int OutputProcessState(const std::string &file_prefix, const std::vector<std::tuple<std::string, double>> &reactions)
{
  fs::path file_path = file_prefix + "process_state.csv";
  std::ofstream csv_file(file_path);
  if (!csv_file.is_open())
  {
    std::cerr << "Failed to open the CSV file " << file_path << "." << std::endl;
    return 1;
  }
  csv_file << "Reaction Label,Rate Constant" << std::endl;
  for (auto &reaction : reactions)
    csv_file << std::get<0>(reaction) << "," << std::setprecision(10)
            << std::get<1>(reaction) << std::endl;
  csv_file.close();
  return 0;
}

// https://www.fluentcpp.com/2020/02/11/reverse-for-loops-in-cpp/
template <typename T>
class reverse
{
private:
  T &iterable_;

public:
  explicit reverse(T &iterable) : iterable_{iterable} {}
  auto begin() const { return std::rbegin(iterable_); }
  auto end() const { return std::rend(iterable_); }
};

// use the reactants as a label for the reaction
// allows us to more easily match reaction information from cam-chem
std::string reactants_to_string(const micm::Process &process)
{
  std::string label{};

  if (auto constant = dynamic_cast<UserDefinedRateConstant *>(process.rate_constant_.get()))
  {
    label = constant->parameters_.label_;
    label += "|";
  }

  // cam-chem seems to use labels that are in reverse
  for (auto &reactant : process.reactants_)
  {
    label += reactant.name_ + "_";
  }
  // delete the trailing underscore
  label.erase(label.size() - 1);

  // Some custom treatments
  if (label == "HCL_O1D" || label == "O1D_HCL")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "CL")
        return "O1D_HCLa";
      if (product.first.name_ == "CLO")
        return "O1D_HCLb";
    }
  if (label == "HBR_O1D" || label == "O1D_HBR")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "BR")
        return "O1D_HBRa";
      if (product.first.name_ == "BRO")
        return "O1D_HBRb";
    }
  if (label == "USER.usr_DMS_OH|DMS_OH" || label == "USER.usr_DMS_OH|OH_DMS")
  {
    if (process.products_.size() == 2)
    {
      return "usr_DMS_OH";
    }
  }
  if (label == "DMS_OH" || label == "OH_DMS")
  {
    if (process.products_.size() == 1)
    {
      return "DMS_OHa";
    }
  }
  if (label == "USER.usr_CO_OH|CO_OH" || label == "USER.usr_CO_OH|OH_CO")
  {
    return "usr_CO_OH";
  }
  if (label == "CLO_CLO")
  {
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "CL2")
        return "CLO_CLOb";
      if (product.first.name_ == "OCLO")
        return "CLO_CLOc";
      if (product.first.name_ == "CL2O2")
        return "tag_CLO_CLO_M";
    }
    return "CLO_CLOa";
  }
  if (label == "BRO_CLO" || label == "BRO_CLO")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "OCLO")
        return "BRO_CLOa";
      if (product.first.name_ == "CL")
        return "BRO_CLOb";
      if (product.first.name_ == "BRCL")
        return "BRO_CLOc";
    }
  if (label == "O1D_N2O" || label == "N2O_O1D")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "NO")
        return "O1D_N2Oa";
      if (product.first.name_ == "N2")
        return "O1D_N2Ob";
    }
  if (label == "CH4_O1D" || label == "O1D_CH4")
  {
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "CH3O2")
        return "O1D_CH4a";
      if (product.first.name_ == "HO2")
        return "O1D_CH4b";
    }
    return "O1D_CH4c";
  }
  if (label == "CH3O2_CH3O2")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "HO2")
        return "CH3O2_CH3O2a";
      if (product.first.name_ == "CH3OH")
        return "CH3O2_CH3O2b";
    }

  if (label == "CL_HO2" || label == "HO2_CL")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "HCL")
        return "CL_HO2a";
      if (product.first.name_ == "CLO")
        return "CL_HO2b";
    }
  if (label == "CLO_OH" || label == "OH_CLO")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "HO2")
        return "CLO_OHa";
      if (product.first.name_ == "HCL")
        return "CLO_OHb";
    }
  if (label == "MTERP_O3" || label == "O3_MTERP")
  {
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "TERPROD1")
        return "MTERP_O3";
    }
    return "MTERP_O3_vbs";
  }
  if (label == "N_NO2" || label == "NO2_N")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "N2O")
        return "N_NO2a";
      if (product.first.name_ == "NO")
        return "N_NO2b";
      if (product.first.name_ == "N2")
        return "N_NO2c";
    }
  if (label == "ENEO2_NO" || label == "NO_ENEO2")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "CH3CHO")
        return "ENEO2_NO";
      if (product.first.name_ == "HONITR")
        return "ENEO2_NOb";
    }
  if (label == "OH_OH")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "H2O")
        return "OH_OH";
      if (product.first.name_ == "H2O2")
        return "OH_OH_M";
    }
  if (label == "O3_O1D" || label == "O1D_O3")
  {
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "O")
        return "O1D_O3a";
    }
    return "O1D_O3";
  }
  if (label == "MTERP_OH" || label == "OH_MTERP")
  {
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "TERPO2")
        return "MTERP_OH";
    }
    return "MTERP_OH_vbs";
  }
  if (label == "MTERP_NO3" || label == "NO3_MTERP")
  {
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "NTERPO2")
        return "MTERP_NO3";
    }
    return "MTERP_NO3_vbs";
  }
  if (label == "ALKO2_NO" || label == "NO_ALKO2")
  {
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "CH3CHO")
        return "ALKO2_NO";
    }
    return "ALKO2_NOb";
  }
  if (label == "BCARY_O3" || label == "O3_BCARY")
  {
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "TERPROD1")
        return "BCARY_O3";
    }
    return "BCARY_O3_vbs";
  }
  if (label == "BCARY_NO3" || label == "NO3_BCARY")
  {
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "NTERPO2")
        return "BCARY_NO3";
    }
    return "BCARY_NO3_vbs";
  }
  if (label == "BCARY_OH" || label == "OH_BCARY")
  {
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "TERPO2")
        return "BCARY_OH";
    }
    return "BCARY_OH_vbs";
  }
  if (label == "BENZENE_OH" || "OH_BENZENE")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "BENZO2VBS")
        return "BENZENE_OH_vbs";
    }
  if (label == "ISOP_OH" || label == "OH_ISOP")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "ISOPO2VBS")
        return "ISOP_OH_vbs";
    }
  if (label == "H_HO2" || label == "HO2_H")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "H2")
        return "H_HO2";
      if (product.first.name_ == "OH")
        return "H_HO2a";
      if (product.first.name_ == "H2O")
        return "H_HO2b";
    }
  if (label == "TOLUENE_OH" || label == "OH_TOLUENE")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "TOLUO2VBS")
        return "TOLUENE_OH_vbs";
    }
  if (label == "NO2_O" || label == "NO2_O")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "NO")
        return "NO2_O";
      if (product.first.name_ == "NO3")
        return "NO2_O_M";
    }
  if (label == "ISOP_O3" || label == "O3_ISOP")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "MVK")
        return "ISOP_O3";
      if (product.first.name_ == "SOAG3")
        return "ISOP_O3_vbs";
    }
  if (label == "ISOP_NO3" || label == "NO3_ISOP")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "ISOPNO3")
        return "ISOP_NO3";
      if (product.first.name_ == "SOAG3")
        return "ISOP_NO3_vbs";
    }
  if (label == "XLYENES_OH" || label == "OH_XYLENES")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "XYLOL")
        return "XYLENES_OH";
      if (product.first.name_ == "XYLEO2VBS")
        return "XYLENES_OH_vbs";
    }
  if (label == "MACRO2_NO" || label == "NO_MACRO2")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "NO2")
        return "MACRO2_NOa";
      if (product.first.name_ == "HONITR")
        return "MACRO2_NOb";
    }
  if (label == "N2O5")
    for (auto &product : process.products_)
    {
      if (product.first.name_ == "NO2")
        return "usr_N2O5_M";
      if (product.first.name_ == "HNO3")
        return "usr_N2O5_aer";
    }
  if (label == "CL2O2")
    for (auto &product : process.products_)
      if (product.first.name_ == "CLO")
          return "usr_CL2O2_M";
  if (label == "HO2_HO2" || label == "HO2_HO2_M" || label == "M_HO2_HO2" ||
      label == "H2O_HO2_HO2" || label == "HO2_HO2_H2O" || label == "H2O_HO2_HO2_M" ||
      label == "M_HO2_HO2_H2O") return "usr_HO2_HO2";

  label += "|";
  for (auto &reactant : reverse(process.reactants_))
  {
    label += reactant.name_ + "_";
  }
  auto pos = label.find(",");
  if (pos != std::string::npos)
  {
    label.erase(pos, 1);
  }

  // delete the trailing underscore
  label.erase(label.size() - 1);

  return label;
}

int OutputDiagnosticData(const std::string file_prefix, const SolverParameters &solver_parameters, InitialConditions initial_conditions)
{
  RosenbrockSolver<> solver{
    solver_parameters.system_,
    solver_parameters.processes_,
    RosenbrockSolverParameters::three_stage_rosenbrock_parameters(1)};

  State state = solver.GetState();
  SetState(initial_conditions, 1, state);
  solver.UpdateState(state);
  Matrix<double> initial_forcing(state.variables_.NumRows(), state.variables_[0].size(), 0.0);
  solver.CalculateForcing(state.rate_constants_, state.variables_, initial_forcing);

  std::vector<std::tuple<std::string, double>> reaction_data;
  std::cout << solver_parameters.processes_.size() << ", " << state.rate_constants_[0].size() << std::endl;
  for (int i = 0; i < state.rate_constants_[0].size(); ++i)
  {
    // convert the rate constant to mixing ratios to mimic what cam-chem does in mo_adjrxt
    double conversion = std::pow(initial_conditions.environments["air_density"], solver_parameters.processes_[i].reactants_.size() - 1);
    double rate_constant = state.rate_constants_[0][i];
    rate_constant *= conversion;


    // for SO3 + 2H2O -> H2SO4 + H20, camchem doesn't account for one of the water concentrations until
    // after the rate constant is calculated, so we need to adjust our to see if we match
    double water_count = 0;
    double so3_count = 0;
    // CAM-Chem multiplies the rate constants by constant species concentrations (M, O2, N2)
    // whose values are concentrations unlike all the other species, which are in mixing ratios
    // (M would be ignored since its mixing ratio would be 1.0 mol mol-1, but M is parameterized in micm, and is thus already included in the rate constant calculation)
    for (auto& reactant : solver_parameters.processes_[i].reactants_)
    {
      if (reactant.name_ == "H2O") ++water_count;
      if (reactant.name_ == "SO3") ++so3_count;
      if (reactant.name_ == "M") rate_constant /= state.conditions_[0].air_density_;
      if (reactant.name_ == "O2") rate_constant *= state.variables_[0][state.variable_map_["O2"]] / state.conditions_[0].air_density_;
      if (reactant.name_ == "N2") rate_constant *= state.variables_[0][state.variable_map_["N2"]] / state.conditions_[0].air_density_;
    }

    // The custom HO2 + HO2 rate constant in CAM-Chem has [H2O] baked in
    int n_HO2 = 0;
    int n_H2O = 0;
    for (auto& reactant : solver_parameters.processes_[i].reactants_)
    {
      if (reactant.name_ == "HO2") ++n_HO2;
      if (reactant.name_ == "H2O") ++n_H2O;
    }
    if (n_HO2 == 2 && n_H2O == 1) rate_constant *= state.variables_[0][state.variable_map_["H2O"]] / state.conditions_[0].air_density_;


    // found the  SO3 + 2H2O -> H2SO4 + H20 reaction, divide by the concentration of water
    if (water_count == 2 && so3_count == 1) {
        rate_constant *= state.variables_[0][state.variable_map_["H2O"]] / state.conditions_[0].air_density_;
    }
#if 0
    if (solver_parameters.processes_[i].reactants_[0].name_ == "HO2" &&
        solver_parameters.processes_[i].reactants_[1].name_ == "NO2")
    {
      auto& rc = dynamic_cast<TroeRateConstant&>(*(solver_parameters.processes_[i].rate_constant_));
      PrintTroe(rc, state.conditions_[0], rate_constant);
    }
#endif
    reaction_data.push_back(std::make_tuple(reactants_to_string(solver_parameters.processes_[i]), rate_constant));
  }

  std::vector<std::tuple<std::string, double, double>> species_data;
  auto species_names = solver_parameters.system_.UniqueNames();
  for (int i = 0; i < species_names.size(); ++i)
    species_data.push_back(std::make_tuple(
      species_names[i],
      state.variables_[0][state.variable_map_[species_names[i]]],
      initial_forcing[0][state.variable_map_[species_names[i]]]));
  for (auto &species : solver_parameters.system_.gas_phase_.species_)
    if (species.IsParameterized())
        species_data.push_back(std::make_tuple(
          species.name_, species.parameterize_(state.conditions_[0]),
          0.0));

    //int result_code = OutputSpeciesState(file_prefix, species_data);
    //if (result_code != 0)
    //  return result_code;

    //result_code = OutputProcessState(file_prefix, reaction_data);
    //if (result_code != 0)
    //  return result_code;

    return 0;
  }

template <template <class> class MatrixType, template <class> class SparseMatrixType>
void solveAndWriteResults(
  const SolverParameters &solver_params,
  InitialConditions &initial_conditions,
  const std::string &matrix_ordering,
  std::ofstream &csvFile,
  int grid_cells,
  const std::string& matrix_ordering_type, 
  int id)
{
  using SolverType = RosenbrockSolver<MatrixType, SparseMatrixType>;

  auto params = RosenbrockSolverParameters::three_stage_rosenbrock_parameters(grid_cells);
  SolverType solver{
    solver_params.system_,
    solver_params.processes_,
    params};

  size_t nspecies = solver_params.system_.StateSize();
  for (auto& species : solver_params.system_.gas_phase_.species_)
    if (species.IsParameterized()) ++nspecies;
  size_t nreactions = solver_params.processes_.size();

  //std::string results_file_name = "micm_" + matrix_ordering + "_" + std::to_string(grid_cells) + "_gridcells_" + std::to_string(nspecies) + "_species_" + std::to_string(nreactions) + "_reactions.csv";
  //fs::path results_file_path = results_file_name;
  //std::ofstream results_file(results_file_path);
  //if (!results_file.is_open())
  //{
  //  std::cerr << "Failed to open the concentration results file." << std::endl;
  //  return;
  //}

  State state = solver.GetState();
  SetState(initial_conditions, grid_cells, state);

  double time = 0;
  typename SolverType::SolverResult result;

  auto timestep = initial_conditions.environments["time_step"];

  auto startTime = std::chrono::high_resolution_clock::now();

  MICM_PROFILE_BEGIN_SESSION("Runtime", "Profile-Runtime-" + matrix_ordering_type + ".json", id);
  while (time < timestep)
  {
    result = solver.Solve(timestep - time, state);
    bool b = result.state_ == SolverState::Converged;
    assert(b);
    state.variables_ = result.result_;
    time += result.final_time_;
  }
  MICM_PROFILE_END_SESSION();

  //auto endTime = std::chrono::high_resolution_clock::now();
  //auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime);

  //csvFile << matrix_ordering << "," << grid_cells << "," << nspecies << "," << nreactions << "," << duration.count() << std::endl;

  //results_file << "Species";
  //for (int i{}; i<grid_cells; ++i) results_file << ",Cell_" << i;
  //results_file << std::endl;

  //results_file << std::setprecision(16);
  //auto map = state.variable_map_;
  //for (const auto &pair : map)
  //{
  //  results_file << pair.first;
  //  for (size_t cell = 0; cell < params.number_of_grid_cells_; ++cell)
  //    results_file << "," << state.variables_[cell][pair.second];
  //  results_file << std::endl;
  //}
  //results_file.close();
}

template <class T>
using SparseMatrixParam = micm::SparseMatrix<T>;
template <class T>
using Vector1MatrixParam = micm::VectorMatrix<T, 1>;
template <class T>
using Vector10MatrixParam = micm::VectorMatrix<T, 10>;
template <class T>
using Vector100MatrixParam = micm::VectorMatrix<T, 100>;
template <class T>
using Vector1000MatrixParam = micm::VectorMatrix<T, 1000>;
template <class T>
using Vector1SparseMatrixParam = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template <class T>
using Vector10SparseMatrixParam = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<10>>;
template <class T>
using Vector100SparseMatrixParam = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<100>>;
template <class T>
using Vector1000SparseMatrixParam = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1000>>;

int main(const int argc, const char *argv[])
{
  if (argc < 3)
  {
    std::cout << "Usage: performance_test_micm_performance <path> <matrix-ordering>" << std::endl;
    std::cout << std::endl;
    std::cout << "  path             Path to the folder holding the configuration data" << std::endl;
    std::cout << "  matrix-ordering  Specify the matrix ordering \"standard\" or \"vector\"" << std::endl;
    return 1;
  }

  int gridSizes[] = {1, 10, 100, 1000};

  //fs::path configDir = fs::path(argv[1]);
  std::string configDir{ argv[1] };
  std::string matrix_ordering{argv[2]};

  SolverConfig solverConfig;
  ConfigParseStatus status = solverConfig.ReadAndParse(fs::path(configDir));

  SolverParameters solver_params = solverConfig.GetSolverParams();

  // add third-body species parameterizaton on air density
  for (auto& species : solver_params.system_.gas_phase_.species_)
    if (species.name_ == "M")
      species.parameterize_ = [](const Conditions& c) { return c.air_density_; };
  for (auto& process : solver_params.processes_)
  {
    for (auto& reactant : process.reactants_)
      if (reactant.name_ == "M")
        reactant.parameterize_ = [](const Conditions& c) { return c.air_density_; };
    for (auto& product : process.products_)
      if (product.first.name_ == "M")
        product.first.parameterize_ = [](const Conditions& c) { return c.air_density_; };
  }

  size_t nspecies = solver_params.system_.StateSize();
  size_t nreactions = solver_params.processes_.size();

  std::cout << "Running MICM performance comparison" << std::endl;
  std::cout << "  config dir: " << configDir << std::endl;
  std::cout << "  number of species: " << nspecies << std::endl;
  std::cout << "  number of reactions: " << nreactions << std::endl;
  std::cout << "  matrix ordering: " << matrix_ordering << std::endl;

  auto initial_conditions = read_initial_conditions(configDir + "/initial_conditions.csv");

  const std::string file_prefix = "micm_" + matrix_ordering + "_" + std::to_string(nspecies) + "_species_" + std::to_string(nreactions) + "_reactions_";

  int result_code = OutputDiagnosticData(file_prefix, solver_params, initial_conditions);
  if (result_code != 0)
    return result_code;

  std::string file_name = file_prefix + "timing_results.csv";
  fs::path csvFilePath = file_name;
  std::ofstream csvFile(csvFilePath);
  //if (!csvFile.is_open())
  //{s
  //  std::cerr << "Failed to open the CSV file " << csvFilePath << "." << std::endl;
  //  return 1;
  //}
  //csvFile << "Matrix,Grid cells,Species,Reactions,Solving Time (nanoseconds)" << std::endl;

  if (matrix_ordering == "standard")
  {
      solveAndWriteResults<Matrix, SparseMatrixParam>(solver_params, initial_conditions, "standard", csvFile, 1, "Sparse-1", 10001);
      solveAndWriteResults<Matrix, SparseMatrixParam>(solver_params, initial_conditions, "standard", csvFile, 10, "Sparse-10", 10010);
      solveAndWriteResults<Matrix, SparseMatrixParam>(solver_params, initial_conditions, "standard", csvFile, 100, "Sparse-100", 10100);
      solveAndWriteResults<Matrix, SparseMatrixParam>(solver_params, initial_conditions, "standard", csvFile, 1000, "Sparse-1000", 11000);
  }
  else if (matrix_ordering == "vector")
  {
    //solveAndWriteResults<Vector1MatrixParam, Vector1SparseMatrixParam>
    //  (solver_params, initial_conditions, "vector", csvFile, 1, "Vector-1", 20001);
    //solveAndWriteResults<Vector10MatrixParam, Vector10SparseMatrixParam>
    //  (solver_params, initial_conditions, "vector", csvFile, 10, "Vector-10", 20010);
    //solveAndWriteResults<Vector100MatrixParam, Vector100SparseMatrixParam>
    //  (solver_params, initial_conditions,"vector", csvFile, 100, "Vector-100", 20100);
    solveAndWriteResults<Vector1000MatrixParam, Vector1000SparseMatrixParam>
      (solver_params, initial_conditions, "vector", csvFile, 1000, "Vector-1000", 21000);
  }
  std::cout << "Done!" << std::endl;

  //csvFile.close();
  return 0;
}