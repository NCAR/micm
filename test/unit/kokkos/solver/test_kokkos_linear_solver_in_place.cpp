#include "../../solver/test_linear_solver_in_place_policy.hpp"

#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/solver/linear_solver.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <functional>

using FloatingPointType = micm::Real;

using Group1VectorMatrix = micm::KokkosDenseMatrix<FloatingPointType, 1>;
using Group2VectorMatrix = micm::KokkosDenseMatrix<FloatingPointType, 2>;
using Group3VectorMatrix = micm::KokkosDenseMatrix<FloatingPointType, 3>;
using Group4VectorMatrix = micm::KokkosDenseMatrix<FloatingPointType, 4>;

using Group1SparseVectorMatrix = micm::KokkosSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::KokkosSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::KokkosSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::KokkosSparseMatrix<FloatingPointType, micm::SparseMatrixVectorOrdering<4>>;

TEST(KokkosLinearSolverInPlace, DenseMatrixVectorOrdering)
{
  TestDenseMatrix<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::LinearSolverInPlace<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  TestDenseMatrix<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      micm::LinearSolverInPlace<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  TestDenseMatrix<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::LinearSolverInPlace<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  TestDenseMatrix<
      Group4VectorMatrix,
      Group4SparseVectorMatrix,
      micm::LinearSolverInPlace<Group4VectorMatrix, Group4SparseVectorMatrix>>();
}

TEST(KokkosLinearSolverInPlace, RandomMatrixVectorOrdering)
{
  TestRandomMatrix<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::LinearSolverInPlace<Group1VectorMatrix, Group1SparseVectorMatrix>>(5);
  TestRandomMatrix<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      micm::LinearSolverInPlace<Group2VectorMatrix, Group2SparseVectorMatrix>>(5);
  TestRandomMatrix<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::LinearSolverInPlace<Group3VectorMatrix, Group3SparseVectorMatrix>>(5);
  TestRandomMatrix<
      Group4VectorMatrix,
      Group4SparseVectorMatrix,
      micm::LinearSolverInPlace<Group4VectorMatrix, Group4SparseVectorMatrix>>(5);
}

TEST(KokkosLinearSolverInPlace, VectorOrderingAgnosticToInitialValue)
{
  micm::Real initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto initial_value : initial_values)
  {
    TestExtremeInitialValue<
        Group1VectorMatrix,
        Group1SparseVectorMatrix,
        micm::LinearSolverInPlace<Group1VectorMatrix, Group1SparseVectorMatrix>>(1, initial_value);
    TestExtremeInitialValue<
        Group2VectorMatrix,
        Group2SparseVectorMatrix,
        micm::LinearSolverInPlace<Group2VectorMatrix, Group2SparseVectorMatrix>>(2, initial_value);
    TestExtremeInitialValue<
        Group3VectorMatrix,
        Group3SparseVectorMatrix,
        micm::LinearSolverInPlace<Group3VectorMatrix, Group3SparseVectorMatrix>>(5, initial_value);
    TestExtremeInitialValue<
        Group4VectorMatrix,
        Group4SparseVectorMatrix,
        micm::LinearSolverInPlace<Group4VectorMatrix, Group4SparseVectorMatrix>>(5, initial_value);
  }
}

TEST(KokkosLinearSolverInPlace, DiagonalMatrixVectorOrdering)
{
  TestDiagonalMatrix<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::LinearSolverInPlace<Group1VectorMatrix, Group1SparseVectorMatrix>>(5);
  TestDiagonalMatrix<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      micm::LinearSolverInPlace<Group2VectorMatrix, Group2SparseVectorMatrix>>(5);
  TestDiagonalMatrix<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::LinearSolverInPlace<Group3VectorMatrix, Group3SparseVectorMatrix>>(5);
  TestDiagonalMatrix<
      Group4VectorMatrix,
      Group4SparseVectorMatrix,
      micm::LinearSolverInPlace<Group4VectorMatrix, Group4SparseVectorMatrix>>(5);
}

TEST(KokkosLinearSolverInPlace, VectorDiagonalMarkowitzReordering)
{
  TestMarkowitzReordering<Group1VectorMatrix, Group1SparseVectorMatrix>();
  TestMarkowitzReordering<Group2VectorMatrix, Group2SparseVectorMatrix>();
  TestMarkowitzReordering<Group3VectorMatrix, Group3SparseVectorMatrix>();
  TestMarkowitzReordering<Group4VectorMatrix, Group4SparseVectorMatrix>();
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
