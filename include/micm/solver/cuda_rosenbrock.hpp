#pragma once 
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include<micm/solver/cuda_rosenbrock.cuh>
#include <micm/util/cuda_param.hpp>
#include <micm/process/process.hpp>
#include <micm/process/process_set.hpp> 
#include <micm.solver/linear_solver.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>
#include <micm/util/jacobian.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>

namespace micm{
class CudaRosenbrockSolver : public RosenbrockSolver{
///@brief Default constructor 
CudaRosenbrockSolver(); 

CudaRosenbrockSolver(const System& system, 
                    const std::vector<Process>& processes, 
                    const RosenbrockSolverParameters& parameters)
: RosenbrockSolver(system, processes, parameters){}; 

CudaRosenbrockSolver(const System& system
                    const std::vector<Process> processes, 
                    const RosenbrockSolverParamters& parameters
                    const std::function<linearSolverPolicy (const SparseMatrixPolicy<double>, double)> create_linear_solver,
                    const std::function<ProcessSetPolicy (const std::vector<Process>& , const std::map<std::string, std::size_t>&)> create_process_set)
: RosenbrockSolver(system, processes, parameters, create_linear_solver, create_process_set){}; 

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy, class ProcesseSetPolicy> 
requires(VectorizableSparse<SparseMatrixPolicy<double>>);
void AlphaMinusJacobian(SparseMatrixPolicy<double>& jacobian, const double& alpha) const
{
    CudaSparseMatrixParam sparseMatrix; 
    sparseMatrix.jacobian_ = jacobian.AsVector(); 
    sparseMatrix.jacobian_size_ = jacobian.AsVector().size(); 
    sparseMatrix.n_grid_ = jacobian.size(); 
    for (auto& element : sparseMatrix.jacobian_)
    {
        element = -element; 
    }
    std::vector<size_t> jacobian_diagonal_element = state_parameters_.jacobian_diagonal_elements_; //just pass this in without storing in another vector 
    AlphaMinusJacobianDriver(sparseMatrix,
                            state_parameters_.jacobian_diagonal_elements_, 
                            alpha);
    
    }; 
}//end CudaRosenbrockSolver
}//end micm 