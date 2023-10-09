#include <micm/configure/solver_config.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

template<class T>
using SparseMatrixPolicy = micm::SparseMatrix<T>;

void print_header()
{
  std::cout << std::setw(5) << "time"
            << "," << std::setw(10) << "O2"
            << "," << std::setw(10) << "O3"
            << "," << std::setw(10) << "O"
            << "," << std::setw(10) << "O1D" << std::endl;
}

template<template<class> class T>
void print_state(double time, micm::State<T>& state)
{
  std::ios oldState(nullptr);
  oldState.copyfmt(std::cout);

  std::cout << std::setw(5) << time << "," << std::flush;

  std::cout << std::scientific << std::setw(10) << std::setprecision(2)
            << state.variables_[0][state.variable_map_["O2"]]
            << "," << std::setw(10) << std::setprecision(2)
            << state.variables_[0][state.variable_map_["O3"]]
            << "," << std::setw(10) << std::setprecision(2)
            << state.variables_[0][state.variable_map_["O"]]
            << "," << std::setw(10) << std::setprecision(2)
            << state.variables_[0][state.variable_map_["O1D"]] << std::endl;

  std::cout.copyfmt(oldState);
}

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
  };

  state.SetConcentrations(solver_params.system_, intial_concentration);

  double time_step = 10;  // s
  int nstep = 20;

  print_header();
  print_state(0, state);

  for (int i = 0; i < nstep; ++i) {

    double elapsed_solve_time = 0;

    while (elapsed_solve_time < time_step) {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
      state.variables_[0] = result.result_.AsVector();
    }
    print_state(time_step * (i + 1), state);
  }

  return 0;
}
