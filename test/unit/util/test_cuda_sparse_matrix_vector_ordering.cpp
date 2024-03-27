#include <gtest/gtest.h>

#include <micm/util/cuda_sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>

#include "cuda_matrix_utils.cuh"
#include "test_sparse_matrix_policy.hpp"

TEST(CudaSparseVectorMatrix, ZeroMatrix)
{
  auto matrix = testZeroMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<3>>();
}

TEST(CudaSparseVectorMatrix, ConstZeroMatrix)
{
  const auto matrix = testConstZeroMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<3>>();
}

TEST(CudaSparseVectorMatrix, CopyAssignmentZeroMatrix)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>::create(2)
                     .with_element(0, 0)
                     .with_element(1, 1)
                     .initial_value(0.0);

  micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>> matrix{ builder };

  matrix.CopyToDevice();
  auto param = matrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  EXPECT_EQ(matrix.AsVector().size(), 2);

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  matrix.CopyToHost();

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseVectorMatrix, CopyAssignmentConstZeroMatrix)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>::create(2)
                     .with_element(0, 0)
                     .with_element(1, 1)
                     .initial_value(0.0);

  const micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>> matrix{ builder };

  auto oneMatrix = matrix;
  oneMatrix.CopyToDevice();

  auto param = oneMatrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  EXPECT_EQ(matrix.AsVector().size(), 2);

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseVectorMatrix, CopyAssignmentDeSynchedHostZeroMatrix)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>::create(2)
                     .with_element(0, 0)
                     .with_element(1, 1)
                     .initial_value(0.0);

  micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>> matrix{ builder };

  matrix.CopyToDevice();

  auto param = matrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  auto oneMatrix = matrix;

  EXPECT_EQ(2, matrix.AsVector().size());
  EXPECT_EQ(2, oneMatrix.AsVector().size());

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  matrix.CopyToHost();

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseVectorMatrix, MoveAssignmentConstZeroMatrix)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>::create(2)
                     .with_element(0, 0)
                     .with_element(1, 1)
                     .initial_value(0.0);

  const micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>> matrix{ builder };

  auto oneMatrix = std::move(matrix);

  EXPECT_EQ(2, oneMatrix.AsVector().size());
  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToDevice();
  auto param = oneMatrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseVectorMatrix, MoveAssignmentDeSyncedHostZeroMatrix)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>>::create(2)
                     .with_element(0, 0)
                     .with_element(1, 1)
                     .initial_value(0.0);

  micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrdering<1>> matrix{ builder };

  EXPECT_EQ(2, matrix.AsVector().size());
  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  matrix.CopyToDevice();
  auto param = matrix.AsDeviceParam();
  micm::cuda::AddOneDriver(param);

  for (auto& elem : matrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  auto oneMatrix = std::move(matrix);

  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 0.0);
  }

  oneMatrix.CopyToHost();

  for (auto& elem : oneMatrix.AsVector())
  {
    EXPECT_EQ(elem, 1.0);
  }
}

TEST(CudaSparseVectorMatrix, SetScalar)
{
  testSetScalar<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<3>>();
}

TEST(CudaSparseVectorMatrix, SingleBlockMatrix)
{
  auto matrix = testSingleBlockMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<4>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 12);
    matrix.AsVector()[elem] = 42;
    EXPECT_EQ(matrix.AsVector()[12], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 8);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[8], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 4);
  EXPECT_EQ(matrix.GroupSize(matrix.FlatBlockSize()), 4 * 4);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(CudaSparseVectorMatrix, ConstSingleBlockMatrix)
{
  auto matrix = testConstSingleBlockMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<2>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 6);
    EXPECT_EQ(matrix.AsVector()[6], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 4);
    EXPECT_EQ(matrix.AsVector()[4], 21);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(matrix.FlatBlockSize()), 2 * 4);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);
}

TEST(CudaSparseVectorMatrix, MultiBlockMatrix)
{
  auto matrix = testMultiBlockMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrdering<2>>();

  {
    std::size_t elem = matrix.VectorIndex(0, 2, 3);
    EXPECT_EQ(elem, 4);
    matrix.AsVector()[elem] = 21;
    EXPECT_EQ(matrix.AsVector()[4], 21);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 2, 1);
    EXPECT_EQ(elem, 10);
    matrix.AsVector()[elem] = 31;
    EXPECT_EQ(matrix.AsVector()[10], 31);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(matrix.GroupSize(matrix.FlatBlockSize()), 2 * 4);
  EXPECT_EQ(matrix.NumberOfGroups(4), 2);
}
