#include "../terminator.hpp"

#include <micm/jit/solver/jit_rosenbrock.hpp>
#include <micm/jit/solver/jit_solver_builder.hpp>
#include <micm/jit/solver/jit_solver_parameters.hpp>
#include <micm/process/process.hpp>
#include <micm/system/system.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>
template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = micm::VectorMatrix<T, 4>;

using Group1SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<4>>;

TEST(JitRosenbrockSolver, Terminator)
{
  auto parameters = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters();
  parameters.relative_tolerance_ = 1.0e-8;
  parameters.max_number_of_steps_ = 100000;
  {
    auto builder = micm::JitSolverBuilder_Doolittle<micm::JitRosenbrockSolverParameters, 1>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 1);
  }
  {
    auto builder = micm::JitSolverBuilder_Doolittle<micm::JitRosenbrockSolverParameters, 2>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 2);
  }
  {
    auto builder = micm::JitSolverBuilder_Doolittle<micm::JitRosenbrockSolverParameters, 3>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 3);
  }
  {
    auto builder = micm::JitSolverBuilder_Doolittle<micm::JitRosenbrockSolverParameters, 4>(parameters).SetIgnoreUnusedSpecies(true);
    TestTerminator(builder, 4);
  }
}