#include "../../process/test_process_set_policy.hpp"

#include <micm/process/process_set.hpp>
#include <micm/kokkos/util/kokkos_dense_matrix.hpp>
#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/types.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <random>

using Group1VectorMatrix = micm::KokkosDenseMatrix<micm::Real, 1>;
using Group2VectorMatrix = micm::KokkosDenseMatrix<micm::Real, 2>;
using Group3VectorMatrix = micm::KokkosDenseMatrix<micm::Real, 3>;
using Group4VectorMatrix = micm::KokkosDenseMatrix<micm::Real, 4>;

using Group1SparseVectorMatrix = micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<1>>;
using Group2SparseVectorMatrix = micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<2>>;
using Group3SparseVectorMatrix = micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<3>>;
using Group4SparseVectorMatrix = micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<4>>;

TEST(KokkosProcessSet, VectorMatrix)
{
  TestProcessSet<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::ProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  TestProcessSet<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      micm::ProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  TestProcessSet<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::ProcessSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  TestProcessSet<
      Group4VectorMatrix,
      Group4SparseVectorMatrix,
      micm::ProcessSet<Group4VectorMatrix, Group4SparseVectorMatrix>>();
}

TEST(KokkosRandomProcessSet, Matrix)
{
  TestRandomSystem<Group1VectorMatrix, Group1SparseVectorMatrix, micm::ProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix>>(
      200, 50, 40);
  TestRandomSystem<Group1VectorMatrix, Group1SparseVectorMatrix, micm::ProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix>>(
      300, 30, 20);
  TestRandomSystem<Group1VectorMatrix, Group1SparseVectorMatrix, micm::ProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix>>(
      400, 100, 80);
}

TEST(KokkosProcessSetAlgebraicVariables, CudaMatrix)
{
  TestAlgebraicMasking<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::ProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
  TestAlgebraicMasking<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      micm::ProcessSet<Group2VectorMatrix, Group2SparseVectorMatrix>>();
  TestAlgebraicMasking<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::ProcessSet<Group3VectorMatrix, Group3SparseVectorMatrix>>();
  TestAlgebraicMasking<
      Group4VectorMatrix,
      Group4SparseVectorMatrix,
      micm::ProcessSet<Group4VectorMatrix, Group4SparseVectorMatrix>>();
}

TEST(KokkosProcessSetFiniteDifferenceJacobian, Matrix)
{
  TestProcessSetFiniteDifferenceJacobian<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::ProcessSet<Group1VectorMatrix, Group1SparseVectorMatrix>>();
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
