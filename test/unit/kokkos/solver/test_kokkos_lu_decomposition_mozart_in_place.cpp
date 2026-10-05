#include "../../solver/test_lu_decomposition_in_place_policy.hpp"

#include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
#include <micm/solver/lu_decomposition.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

using Sparse1 = micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<1>>;
using Sparse2 = micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<2>>;
using Sparse3 = micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<3>>;
using Sparse4 = micm::KokkosSparseMatrix<micm::Real, micm::SparseMatrixVectorOrdering<4>>;

using LU1 = micm::LuDecompositionMozartInPlace<Sparse1>;
using LU2 = micm::LuDecompositionMozartInPlace<Sparse2>;
using LU3 = micm::LuDecompositionMozartInPlace<Sparse3>;
using LU4 = micm::LuDecompositionMozartInPlace<Sparse4>;

TEST(KokkosLuDecompositionMozartInPlace, DenseMatrixVectorOrdering)
{
  TestDenseMatrix<Sparse1, LU1>();
  TestDenseMatrix<Sparse2, LU2>();
  TestDenseMatrix<Sparse3, LU3>();
  TestDenseMatrix<Sparse4, LU4>();
}

TEST(KokkosLuDecompositionMozartInPlace, RandomMatrixVectorOrdering)
{
  TestRandomMatrix<Sparse1, LU1>(5);
  TestRandomMatrix<Sparse2, LU2>(5);
  TestRandomMatrix<Sparse3, LU3>(5);
  TestRandomMatrix<Sparse4, LU4>(5);
}

TEST(KokkosLuDecompositionMozartInPlace, DiagonalMatrixVectorOrdering)
{
  TestDiagonalMatrix<Sparse1, LU1>(5);
  TestDiagonalMatrix<Sparse2, LU2>(5);
  TestDiagonalMatrix<Sparse3, LU3>(5);
  TestDiagonalMatrix<Sparse4, LU4>(5);
}

TEST(KokkosLuDecompositionMozartInPlace, AgnosticToInitialValueVectorOrdering)
{
  micm::Real initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };
  for (auto& value : initial_values)
  {
    TestExtremeValueInitialization<Sparse1, LU1>(5, value);
    TestExtremeValueInitialization<Sparse2, LU2>(5, value);
    TestExtremeValueInitialization<Sparse3, LU3>(5, value);
    TestExtremeValueInitialization<Sparse4, LU4>(5, value);
  }
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}