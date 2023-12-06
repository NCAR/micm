#pragma once 

#include <micm/util/cuda_param.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/solver/cuda_rosenbrock.cuh>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <micm/process/process.hpp>
#include <micm/process/process_set.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/solver/rosenbrock_solver_parameters.hpp>
#include <micm/solver/state.hpp>
#include <micm/system/system.hpp>
#include <micm/util/jacobian.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <string>
#include <vector>
#include <micm/process/cuda_process_set.hpp>
#include <micm/solver/cuda_linear_solver.hpp>

namespace micm{

 template<
      template<class> class MatrixPolicy = Matrix,
      template<class> class SparseMatrixPolicy = StandardSparseMatrix,
      class LinearSolverPolicy = CudaLinearSolver<double, SparseMatrixPolicy>,
      class ProcessSetPolicy = CudaProcessSet>

class CudaRosenbrockSolver : public RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>{
///@brief Default constructor 
public:
CudaRosenbrockSolver(); 


CudaRosenbrockSolver(const System& system, 
                    const std::vector<Process>& processes, 
                    const RosenbrockSolverParameters& parameters)
: RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>(system, processes, parameters){}; 


CudaRosenbrockSolver(const System& system,
                    const std::vector<Process> processes, 
                    const RosenbrockSolverParameters& parameters,
                    const std::function<LinearSolverPolicy (const SparseMatrixPolicy<double>, double)> create_linear_solver,
                    const std::function<ProcessSetPolicy (const std::vector<Process>& , const std::map<std::string, std::size_t>&)> create_process_set)
: RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy, ProcessSetPolicy>(system, processes, parameters, create_linear_solver, create_process_set){}; 



void AlphaMinusJacobian(SparseMatrixPolicy<double>& jacobian, double alpha) const
requires VectorizableSparse<SparseMatrixPolicy<double>>
{
    
    for (auto& element : jacobian.AsVector())
    {
        element = -element; 
    }

     CudaSparseMatrixParam sparseMatrix; 
    sparseMatrix.jacobian_ = jacobian.AsVector().data(); 
    sparseMatrix.jacobian_size_ = jacobian.AsVector().size(); 
    sparseMatrix.n_grids_ = jacobian.size(); 
   
    micm::cuda::AlphaMinusJacobianDriver(sparseMatrix,
                            this->state_parameters_.jacobian_diagonal_elements_, 
                            alpha);
    
        }
    }; //end CudaRosenbrockSolver
}//end micm 