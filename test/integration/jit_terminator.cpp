#include "terminator.hpp"

#include <micm/process/process.hpp>
#include <micm/solver/jit_rosenbrock.hpp>
#include <micm/system/system.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

template<std::size_t number_of_grid_cells, template<class> class MatrixPolicy, class SparseMatrixPolicy>
void RunTerminatorTest()
{
  auto solver_params = micm::RosenbrockSolverParameters::ThreeStageRosenbrockParameters(number_of_grid_cells, true);
  solver_params.relative_tolerance_ = 1.0e-8;
  solver_params.max_number_of_steps_ = 100000;

  TestTerminator<
      micm::JitRosenbrockSolver<
          MatrixPolicy,
          SparseMatrixPolicy,
          micm::JitLinearSolver<number_of_grid_cells, SparseMatrixPolicy>,
          micm::JitProcessSet<number_of_grid_cells>>>(
      number_of_grid_cells, solver_params);
}

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
  RunTerminatorTest<1, Group1VectorMatrix, Group1SparseVectorMatrix>();
  RunTerminatorTest<2, Group2VectorMatrix, Group2SparseVectorMatrix>();
  RunTerminatorTest<3, Group3VectorMatrix, Group3SparseVectorMatrix>();
  RunTerminatorTest<4, Group4VectorMatrix, Group4SparseVectorMatrix>();
}