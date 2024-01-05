#include <gtest/gtest.h>

#include <micm/process/process.hpp>
#include <micm/solver/jit_rosenbrock.hpp>
#include <micm/system/system.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include "terminator.hpp"

template<std::size_t number_of_grid_cells, template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy>
void RunTerminatorTest()
{
  TestTerminator<
      MatrixPolicy,
      micm::JitRosenbrockSolver<
          MatrixPolicy,
          SparseMatrixPolicy,
          micm::JitLinearSolver<number_of_grid_cells, SparseMatrixPolicy>,
          micm::JitProcessSet<number_of_grid_cells>>>(
      [&](const micm::System& s, const std::vector<micm::Process>& p)
          -> micm::JitRosenbrockSolver<
              MatrixPolicy,
              SparseMatrixPolicy,
              micm::JitLinearSolver<number_of_grid_cells, SparseMatrixPolicy>,
              micm::JitProcessSet<number_of_grid_cells>>
      {
        auto jit{ micm::JitCompiler::create() };
        if (auto err = jit.takeError())
        {
          llvm::logAllUnhandledErrors(std::move(err), llvm::errs(), "[JIT Error]");
          EXPECT_TRUE(false);
        }
        auto solver_params = micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters(number_of_grid_cells, true);
        solver_params.absolute_tolerance_ = 1.0e-20;
        solver_params.relative_tolerance_ = 1.0e-8;
        solver_params.max_number_of_steps_ = 100000;
        return micm::JitRosenbrockSolver<
            MatrixPolicy,
            SparseMatrixPolicy,
            micm::JitLinearSolver<number_of_grid_cells, SparseMatrixPolicy>,
            micm::JitProcessSet<number_of_grid_cells>>{ jit.get(), s, p, solver_params };
      },
      number_of_grid_cells);
}

template<class T>
using Group1VectorMatrix = micm::VectorMatrix<T, 1>;
template<class T>
using Group2VectorMatrix = micm::VectorMatrix<T, 2>;
template<class T>
using Group3VectorMatrix = micm::VectorMatrix<T, 3>;
template<class T>
using Group4VectorMatrix = micm::VectorMatrix<T, 4>;

template<class T>
using Group1SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<1>>;
template<class T>
using Group2SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<2>>;
template<class T>
using Group3SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<3>>;
template<class T>
using Group4SparseVectorMatrix = micm::SparseMatrix<T, micm::SparseMatrixVectorOrdering<4>>;

TEST(JitRosenbrockSolver, Terminator)
{
  RunTerminatorTest<1, Group1VectorMatrix, Group1SparseVectorMatrix>();
  RunTerminatorTest<2, Group2VectorMatrix, Group2SparseVectorMatrix>();
  RunTerminatorTest<3, Group3VectorMatrix, Group3SparseVectorMatrix>();
  RunTerminatorTest<4, Group4VectorMatrix, Group4SparseVectorMatrix>();
}