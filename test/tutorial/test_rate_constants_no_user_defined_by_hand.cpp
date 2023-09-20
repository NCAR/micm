#include <iomanip>
#include <iostream>
#include <map>

// Each rate constant is in its own header file
#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/process/branched_rate_constant.hpp>
#include <micm/process/surface_rate_constant.hpp>
#include <micm/process/ternary_chemical_activation_rate_constant.hpp>
#include <micm/process/troe_rate_constant.hpp>
#include <micm/process/tunneling_rate_constant.hpp>
#include <micm/solver/rosenbrock.hpp>

// Use our namespace so that this example is easier to read
using namespace micm;

// The Rosenbrock solver can use many matrix ordering types
// Here, we use the default ordering, but we still need to provide a templated
// Arguent to the solver so it can use the proper ordering with any data type
template<class T>
using SparseMatrixPolicy = SparseMatrix<T>;


void print_header(){
  std::cout 
      << std::setw(8) << "time" << ", "
      << std::setw(22) << "A" << ", "
      << std::setw(22) << "B" << ", "
      << std::setw(22) << "C" << ", "
      << std::setw(22) << "D" 
      << std::endl;
}

template<template<class> class T>
void print_state(double time, State<T>& state) {
  std::cout 
          << std::defaultfloat
          << std::setw(8) << time << ", "
          << std::setw(20) << std::setprecision(16) 
          << std::setw(20) << std::scientific 
          << std::setw(20) << state.variables_[0][state.variable_map_["A"]] << ", "
          << std::setw(20) << state.variables_[0][state.variable_map_["B"]] << ", " 
          << std::setw(20) << state.variables_[0][state.variable_map_["C"]] << ", " 
          << std::setw(20) << state.variables_[0][state.variable_map_["D"]] 
          << std::endl;
}

int main(const int argc, const char *argv[])
{
  auto a = Species("A");
  auto b = Species("B");
  auto c = Species(
      "C",
      std::map<std::string, double>{ { "molecular weight [kg mol-1]", 0.025 },
                                     { "diffusion coefficient [m2 s-1]", 2.3e2 } });
  auto d = Species("D");

  Phase gas_phase{ std::vector<Species>{ a, b, c, d } };

  Process r1 = Process::create()
                   .reactants({ a })
                   .products({ yields(b, 1) })
                   .rate_constant(ArrheniusRateConstant({ .A_ = 2.15e-11, .B_ = 0, .C_ = 110 }))
                   .phase(gas_phase);

  // a branched reaction has two output pathways
  // this is represnted internal to micm as two different reactions
  auto branched_params = BranchedRateConstantParameters{ .X_ = 1.2, .Y_ = 204.3, .a0_ = 1.0e-3, .n_ = 2 };
  branched_params.branch_ = BranchedRateConstantParameters::Branch::Alkoxy;

  Process r2 = Process::create()
                   .reactants({ b })
                   .products({ yields(c, 1) })
                   .rate_constant(BranchedRateConstant(branched_params))
                   .phase(gas_phase);

  branched_params.branch_ = BranchedRateConstantParameters::Branch::Nitrate;
  Process r3 = Process::create()
                   .reactants({ b })
                   .products({ yields(d, 1) })
                   .rate_constant(BranchedRateConstant(branched_params))
                   .phase(gas_phase);

  // to have a stoichiemetric coefficient of more than one for reactants,
  // list the reactant that many times
  // A surface rate constant also needs to know the effective radius and particle number concentration
  // we will set those later
  Process r4 = Process::create()
                   .reactants({ c, c })
                   .products({ yields(a, 1) })
                   .rate_constant(SurfaceRateConstant({ .label_ = "c", .species_ = c, .reaction_probability_ = 0.74 }))
                   .phase(gas_phase);

  Process r5 = Process::create()
                   .reactants({ d })
                   .products({ yields(b, 2) })
                   .rate_constant(TernaryChemicalActivationRateConstant({ .k0_A_ = 1.2,
                                                                          .k0_B_ = 2.3,
                                                                          .k0_C_ = 302.3,
                                                                          .kinf_A_ = 2.6,
                                                                          .kinf_B_ = -3.1,
                                                                          .kinf_C_ = 402.1,
                                                                          .Fc_ = 0.9,
                                                                          .N_ = 1.2 }))
                   .phase(gas_phase);

  Process r6 = Process::create()
                   .reactants({ d })
                   .products({ yields(b, 2) })
                   .rate_constant(TunnelingRateConstant({ .A_ = 1.2, .B_ = 2.3, .C_ = 302.3 }))
                   .phase(gas_phase);

  auto chemical_system = System(micm::SystemParameters{ .gas_phase_ = gas_phase });
  auto reactions = std::vector<micm::Process>{ r1, r2, r3, r4, r5, r6 };

  RosenbrockSolver<Matrix, SparseMatrixPolicy> solver{ chemical_system,
                                                       reactions,
                                                       RosenbrockSolverParameters::three_stage_rosenbrock_parameters() };
  State state = solver.GetState();

  std::cout << r4.rate_constant_->CustomParameters()[0] << std::endl;

  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa

  state.SetConcentration(a, 1.0);  // mol m-3
  state.SetConcentration(b, 0.0);  // mol m-3
  state.SetConcentration(c, 0.0);  // mol m-3
  state.SetConcentration(d, 0.0);  // mol m-3

  state.SetCustomRateParameter("c.effective radius [m]", 1e-7);
  state.SetCustomRateParameter("c.particle number concentration [# m-3]", 2.5e6);

  // choose a timestep a print the initial state
  double time_step = 500;

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
      // std::cout << "solver state: " << StateToString(result.state_) << std::endl;
      state.variables_[0] = result.result_.AsVector();
    }

    print_state(time_step*(i+1), state);
  }

  return 0;
}
