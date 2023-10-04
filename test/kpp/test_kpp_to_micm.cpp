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

  // from Seinfeld and Pandas 3rd ed. table 5.1 p.121, z = 30 km
  state.conditions_[0].temperature_ = 227.0;  // K
  state.conditions_[0].pressure_ = 1200.0;    // Pa

  std::unordered_map<std::string, std::vector<double>> intial_concentration = {
    { "M",   { 3.1e17 } },  // molecules cm-3, S&P3e table 5.1, z = 30 km
    { "O2",  { 6.5e16 } },  // [O2] ~ 0.21 [M]
    { "O3",  { 2.0e12 } },  // typical [O3] mid-latitude z ~ 30 km
    { "O",   { 6.0e7 } },   // [O] / [O3] ~ 3e-5, S&P3e p.124
    { "O1D", { 0.0 } },     //
    { "NO2", { 8.0e8 } },   // ~ 8 ppb
    { "NO",  { 4.0e8 } },   // ~ 4 ppb
  };

  state.SetConcentrations(solver_params.system_, intial_concentration);

  double time_step = 60;  // s
  int nstep = 3;

  for (int i = 0; i < nstep; ++i) {

    double elapsed_solve_time = 0;

    std::cout << "iteration, elapsed_time: "
        << i << ", " << elapsed_solve_time << std::endl;

    while (elapsed_solve_time < time_step) {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
      state.variables_[0] = result.result_.AsVector();
    }
  }

  return 0;
}
