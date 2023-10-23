#include <micm/configure/solver_config.hpp>
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/system/species.hpp>

template<class T>
using SparseMatrixPolicy = micm::SparseMatrix<T>;

void print_header()
{
  std::cout << std::setw(5) << "time"
            << "," << std::setw(10) << "M"
            << "," << std::setw(11) << "O2"
            << "," << std::setw(11) << "O3"
            << "," << std::setw(11) << "O"
            << "," << std::setw(11) << "O1D" << std::endl;
}

template<template<class> class T, template<class> class D>
void print_state(double time, micm::State<T, D>& state)
{
  std::ios oldState(nullptr);
  oldState.copyfmt(std::cout);

  std::cout << std::setw(5) << time << "," << std::flush;

  std::cout << std::scientific << std::setw(10) << std::setprecision(3) << state.variables_[0][state.variable_map_["M"]]
            << "," << std::setw(11) << std::setprecision(3) << state.variables_[0][state.variable_map_["O2"]] << ","
            << std::setw(11) << std::setprecision(3) << state.variables_[0][state.variable_map_["O3"]] << ","
            << std::setw(11) << std::setprecision(3) << state.variables_[0][state.variable_map_["O"]] << "," << std::setw(11)
            << std::setprecision(3) << state.variables_[0][state.variable_map_["O1D"]] << std::endl;

  std::cout.copyfmt(oldState);
}

int main(const int argc, const char* argv[])
{
  micm::SolverConfig solver_config;

  // Read and parse the configure files
  micm::ConfigParseStatus status = solver_config.ReadAndParse("./configs/kpp_chapman");

  if (status != micm::ConfigParseStatus::Success)
  {
    throw "Parsing failed";
  }

  micm::SolverParameters solver_params = solver_config.GetSolverParams();

  auto& process_vector = solver_params.processes_;

  // Print reactions from reactions configuration
  for (int i = 0; i < process_vector.size(); i++)
  {
    int n_reactants = process_vector[i].reactants_.size();
    for (int j = 0; j < n_reactants - 1; j++)
    {
      std::cout << process_vector[i].reactants_[j].name_ << " + ";
    }
    std::cout << process_vector[i].reactants_[n_reactants - 1].name_;

    std::cout << " --> ";

    int n_products = process_vector[i].products_.size();
    for (int j = 0; j < n_products - 1; j++)
    {
      std::cout << process_vector[i].products_[j].first.name_ << " + ";
    }
    std::cout << process_vector[i].products_[n_products - 1].first.name_;

    std::vector<std::string> param_labels = process_vector[i].rate_constant_->CustomParameters();
    for (int k = 0; k < param_labels.size(); k++)
    {
      std::cout << " " << param_labels[k];
    }

    std::cout << std::endl;
  }

  auto chemical_system = solver_params.system_;
  auto reactions = solver_params.processes_;

  micm::RosenbrockSolver<micm::Matrix, SparseMatrixPolicy> solver{
    chemical_system, reactions, micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters()
  };

  micm::State state = solver.GetState();

  state.conditions_[0].temperature_ = 227.0;  // K
  state.conditions_[0].pressure_ = 1200.0;    // Pa

  state.SetCustomRateParameter("PHOTO.R1", 6.0e-11);  // s^-1 j_O2
  state.SetCustomRateParameter("PHOTO.R3", 1.0e-3);   // s^-1 j_O3
  state.SetCustomRateParameter("PHOTO.R5", 1.0e-3);   // s^-1 j_O3

  // Define initial concentrations from Seinfeld & Pandis 3e
  //   and convert to SI units
  double N_Avogadro = 6.02214076e23;
  // molecules cm-3 -> mol m-3, z = 30 km, S&P3e table 5.1 p. 121
  double n_M = 3.1e17 * 1.0e6 / N_Avogadro;
  double n_O2 = 0.21 * n_M;  // [O2] ~ 0.21 [M]
  // typical [O3] mid-latitude z ~ 30 km
  double n_O3 = 2.0e12 * 1.0e6 / N_Avogadro;
  // [O] / [O3] ~ 3e-5, S&P3e p. 124

  micm::Species M("M");
  micm::Species O2("O2");
  micm::Species O3("O3");
  micm::Species O("O");
  micm::Species O1D("O1D");

  state.SetConcentration(M, n_M);
  state.SetConcentration(O2, n_O2);
  state.SetConcentration(O3, n_O3);
  state.SetConcentration(O, 0.0);
  state.SetConcentration(O1D, 0.0);

  double time_step = 3600;  // s
  int nstep = 24;

  print_header();
  print_state(0, state);

  for (int i = 0; i < nstep; ++i)
  {
    double elapsed_solve_time = 0;

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
      state.variables_[0] = result.result_.AsVector();
    }
    print_state(time_step * (i + 1), state);
  }

  return 0;
}
