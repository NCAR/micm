#include <micm/configure/solver_config.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

template<class T>
using SparseMatrixPolicy = micm::SparseMatrix<T>;

int main(const int argc, const char *argv[])
{

  micm::SolverConfig solver_config;

  // Read and parse the configure files
  micm::ConfigParseStatus status = solver_config.ReadAndParse(
    "./configs/kpp_chapman");

  if (status != micm::ConfigParseStatus::Success)
  {
    throw "Parsing failed";
  }

  micm::SolverParameters solver_params = solver_config.GetSolverParams();

  auto& process_vector = solver_params.processes_;

  for (int i = 0; i < process_vector.size(); i++) {

    int n_reactants = process_vector[i].reactants_.size();
    for (int j = 0; j < n_reactants - 1; j++) {
      std::cout << process_vector[i].reactants_[j].name_ << " + ";
    }
    std::cout << process_vector[i].reactants_[n_reactants - 1].name_;

    std::cout << " --> ";

    int n_products = process_vector[i].products_.size();
    for (int j = 0; j < n_products - 1; j++) {
      std::cout << process_vector[i].products_[j].first.name_ << " + ";
    }
    std::cout << process_vector[i].products_[n_products - 1].first.name_;

    std::cout << std::endl;
  }

  auto chemical_system = solver_params.system_;
  auto reactions = solver_params.processes_;

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixPolicy> solver{
    chemical_system, reactions,
    micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };

  micm::State state = solver.GetState();

  state.conditions_[0].temperature_ = 217.0;  // K
  state.conditions_[0].pressure_ = 30000.0;   // Pa

  return 0;
}
