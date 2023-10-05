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

  double unit_conversion = 1.0e6 / 6.023e23; // convert molecules cm-3 to mol m-3

  // from Seinfeld and Pandas 3rd ed. table 5.1 p.121, z = 30 km
  state.conditions_[0].temperature_ = 227.0;  // K
  state.conditions_[0].pressure_ = 1200.0;    // Pa

  double N_Avogadro = 6.02214076e23;
  // molecules cm-3 -> mol m-3, S&P3e table 5.1, z = 30 km
  double n_M = 3.1e17 * 1.0e6 / N_Avogadro;
  // typical [O3] mid-latitude z ~ 30 km
  double n_O3 = 2.0e12 * 1.0e6 / N_Avogadro;

  std::unordered_map<std::string, std::vector<double>> intial_concentration = {
    { "M",   { n_M} },
    { "O2",  { 0.21 * n_M } },  // [O2] ~ 0.21 [M]
    { "O3",  { n_O3 } },
    { "O",   { 3.0e-5 * n_O3 } },  // [O] / [O3] ~ 3e-5, S&P3e p.124
    { "O1D", { 0.0 } },
    { "NO2", { 8.0e-9 * n_M } },  // ~ 8 ppb
    { "NO",  { 4.0e-9 * n_M} },   // ~ 4 ppb
  };

  state.SetConcentrations(solver_params.system_, intial_concentration);

  double time_step = 600;  // s
  int nstep = 20;

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
