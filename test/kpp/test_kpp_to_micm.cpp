#include <micm/configure/solver_config.hpp>

int main(const int argc, const char *argv[])
{

  micm::SolverConfig solver_config;

  // Read and parse the configure files
  micm::ConfigParseStatus status = solver_config.ReadAndParse(
    "etc/configs/micm");

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

  return 0;
}
