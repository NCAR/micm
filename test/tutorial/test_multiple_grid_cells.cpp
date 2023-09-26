#include <iomanip>
#include <iostream>
#include <map>
#include <algorithm>

// Each rate constant is in its own header file
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/process/user_defined_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

// Use our namespace so that this example is easier to read
using namespace micm;

// The Rosenbrock solver can use many matrix ordering types
// Here, we use the default ordering, but we still need to provide a templated
// Arguent to the solver so it can use the proper ordering with any data type
template<class T>
using SparseMatrixPolicy = SparseMatrix<T>;

void print_header()
{
  std::cout << std::setw(5) << "time"
            << "," << std::setw(5) << "grid"
            << "," << std::setw(10) << "O"
            << "," << std::setw(10) << "O1D"
            << "," << std::setw(10) << "O2"
            << "," << std::setw(10) << "O3"
            << "," << std::setw(10) << "M"
            << "," << std::setw(10) << "Ar"
            << "," << std::setw(10) << "N2" 
            << "," << std::setw(10) << "H2O" 
            << "," << std::setw(10) << "CO2" 
            << std::endl;
}

template<template<class> class T>
void print_state(double time, State<T>& state)
{
  std::ios oldState(nullptr);
  oldState.copyfmt(std::cout);


  std::cout << std::setw(5) << time << ",";
  std::cout << std::scientific << std::setprecision(2)
            << std::setw(6) << "1,"
            << std::setw(10) << state.variables_[0][state.variable_map_["O"]] << "," 
            << std::setw(10) << state.variables_[0][state.variable_map_["O1D"]] << ","
            << std::setw(10) << state.variables_[0][state.variable_map_["O2"]] << "," 
            << std::setw(10) << state.variables_[0][state.variable_map_["O3"]] << "," 
            << std::setw(10) << state.variables_[0][state.variable_map_["M"]] << "," 
            << std::setw(10) << state.variables_[0][state.variable_map_["Ar"]] << "," 
            << std::setw(10) << state.variables_[0][state.variable_map_["N2"]] << "," 
            << std::setw(10) << state.variables_[0][state.variable_map_["H2O"]] << "," 
            << std::setw(10) << state.variables_[0][state.variable_map_["CO2"]]
            << std::endl;

  std::cout.copyfmt(oldState);
  std::cout << std::setw(5) << time << ",";
  std::cout << std::scientific << std::setprecision(2)
            << std::setw(6) << "2,"
            << std::setw(10) << state.variables_[1][state.variable_map_["O"]] << "," 
            << std::setw(10) << state.variables_[1][state.variable_map_["O1D"]] << ","
            << std::setw(10) << state.variables_[1][state.variable_map_["O2"]] << "," 
            << std::setw(10) << state.variables_[1][state.variable_map_["O3"]] << "," 
            << std::setw(10) << state.variables_[1][state.variable_map_["M"]] << "," 
            << std::setw(10) << state.variables_[1][state.variable_map_["Ar"]] << "," 
            << std::setw(10) << state.variables_[1][state.variable_map_["N2"]] << "," 
            << std::setw(10) << state.variables_[1][state.variable_map_["H2O"]] << "," 
            << std::setw(10) << state.variables_[1][state.variable_map_["CO2"]]
            << std::endl;

  std::cout.copyfmt(oldState);
  std::cout << std::setw(5) << time << ",";
  std::cout << std::scientific << std::setprecision(2)
            << std::setw(6) << "3,"
            << std::setw(10) << state.variables_[2][state.variable_map_["O"]] << "," 
            << std::setw(10) << state.variables_[2][state.variable_map_["O1D"]] << ","
            << std::setw(10) << state.variables_[2][state.variable_map_["O2"]] << "," 
            << std::setw(10) << state.variables_[2][state.variable_map_["O3"]] << "," 
            << std::setw(10) << state.variables_[2][state.variable_map_["M"]] << "," 
            << std::setw(10) << state.variables_[2][state.variable_map_["Ar"]] << "," 
            << std::setw(10) << state.variables_[2][state.variable_map_["N2"]] << "," 
            << std::setw(10) << state.variables_[2][state.variable_map_["H2O"]] << "," 
            << std::setw(10) << state.variables_[2][state.variable_map_["CO2"]]
            << std::endl;

  std::cout.copyfmt(oldState);
}

int main(){
  auto o =   Species("O");
  auto o1d = Species("O1D");
  auto o2 =  Species("O2");
  auto o3 =  Species("O3");
  auto m =   Species("M");
  auto ar =  Species("Ar");
  auto n2 =  Species("N2");
  auto h2o = Species("H2O");
  auto co2 = Species("CO2");

  Phase gas_phase{ std::vector<Species>{ o, o1d, o2, o3, m, ar, n2, h2o, co2 } };

  Process r1 = Process::create()
                         .reactants({ o1d, n2 })
                         .products({ yields(o, 1), yields(n2, 1) })
                         .rate_constant(ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                         .phase(gas_phase);

  Process r2 = Process::create()
                         .reactants({ o1d, o2 })
                         .products({ yields(o, 1), yields(o2, 1) })
                         .rate_constant(ArrheniusRateConstant(
                             ArrheniusRateConstantParameters{ .A_ = 3.3e-11, .B_ = 0, .C_ = 55 }))
                         .phase(gas_phase);

  Process r3 = Process::create()
                         .reactants({ o, o3 })
                         .products({ yields(o2, 2) })
                         .rate_constant(ArrheniusRateConstant(
                             ArrheniusRateConstantParameters{ .A_ = 8e-12, .B_ = 0, .C_ = -2060 }))
                         .phase(gas_phase);

  Process r4 = Process::create()
                         .reactants({ o, o2, m })
                         .products({ yields(o3, 1), yields(m, 1) })
                         .rate_constant(ArrheniusRateConstant(
                             ArrheniusRateConstantParameters{ .A_ = 6.0e-34, .B_ = 0, .C_ = 2.4 }))
                         .phase(gas_phase);

  Process photo_1 = Process::create()
                              .reactants({ o2 })
                              .products({ yields(o, 2) })
                              .rate_constant(UserDefinedRateConstant({ .label_ = "jO2" }))
                              .phase(gas_phase);

  Process photo_2 = Process::create()
                              .reactants({ o3 })
                              .products({ yields(o1d, 1), yields(o2, 1) })
                              .rate_constant(UserDefinedRateConstant({ .label_ = "jO3a" }))
                              .phase(gas_phase);

  Process photo_3 = Process::create()
                              .reactants({ o3 })
                              .products({ yields(o, 1), yields(o2, 1) })
                              .rate_constant(UserDefinedRateConstant({ .label_ = "jO3b" }))
                              .phase(gas_phase);

  RosenbrockSolver<Matrix, SparseMatrixPolicy> solver{
    System(SystemParameters{ .gas_phase_ = gas_phase }),
    std::vector<Process>{ r1, r2, r3, r4, photo_1, photo_2, photo_3 },
    RosenbrockSolverParameters::three_stage_rosenbrock_parameters(3, false)
  };

  State<Matrix> state = solver.GetState();

  std::vector<double> concentrations{ 1, 1, 1, 2, 2, 2, 3, 3, 3 };
  state.variables_[0] = concentrations;

  std::transform(concentrations.begin(), concentrations.end(), concentrations.begin(), [](auto& c){return c*2;});
  state.variables_[1] = concentrations;

  std::transform(concentrations.begin(), concentrations.end(), concentrations.begin(), [](auto& c){return c/4;});
  state.variables_[2] = concentrations;

  state.SetCustomRateParameter("jO2",  {0.01, 0.01, 0.01});
  state.SetCustomRateParameter("jO3a", {0.02, 0.02, 0.02});
  state.SetCustomRateParameter("jO3b", {0.03, 0.03, 0.03});

  state.conditions_[0].temperature_ = 284.19;  // [K]
  state.conditions_[0].pressure_ = 101245.0;   // [Pa]
  state.conditions_[1].temperature_ = 215.02;  // [K]
  state.conditions_[1].pressure_ = 100789.2;   // [Pa]
  state.conditions_[2].temperature_ = 299.31;  // [K]
  state.conditions_[2].pressure_ = 101398.0;   // [Pa]

  // choose a timestep and print the initial state
  double time_step = 500;  // s

  print_header();
  print_state(0, state);

  // solve for ten iterations
  for (int i = 0; i < 10; ++i)
  {
    // Depending on how stiff the system is
    // the solver integration step may not be able to solve for the full time step
    // so we need to track how much time the solver was able to integrate for and continue
    // solving until we finish
    double elapsed_solve_time = 0;

    while (elapsed_solve_time < time_step)
    {
      auto result = solver.Solve(time_step - elapsed_solve_time, state);
      elapsed_solve_time = result.final_time_;
      state.variables_[0] = result.result_.AsVector();
    }

    print_state(time_step * (i + 1), state);
  }
}