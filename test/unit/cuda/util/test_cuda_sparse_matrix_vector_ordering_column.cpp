#include "../../util/test_sparse_matrix_policy.hpp"
#include "cuda_matrix_utils.cuh"

#include <micm/cuda/util/cuda_sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering_compressed_sparse_column.hpp>
#include <gtest/gtest.h>

/* Below are the policy tests on the CPU */

TEST(SparseVectorCompressedColumnMatrix, ZeroMatrix)
{
  testZeroMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
}

TEST(SparseVectorCompressedColumnMatrix, ConstZeroMatrix)
{
  testConstZeroMatrix<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorCompressedColumnMatrix, SetScalar)
{
  testSetScalar<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
}

TEST(SparseVectorCompressedColumnMatrix, AddToDiagonal)
{
  testAddToDiagonal<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testAddToDiagonal<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testAddToDiagonal<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, Print)
{
  testPrint<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testPrint<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testPrint<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testPrint<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

TEST(SparseVectorCompressedColumnMatrix, PrintNonZero)
{
  testPrintNonZero<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>();
  testPrintNonZero<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<2>>();
  testPrintNonZero<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<3>>();
  testPrintNonZero<micm::CudaSparseMatrix, micm::SparseMatrixVectorOrderingCompressedSparseColumn<4>>();
}

/* These are the customized tests for GPU only */

TEST(CudaSparseMatrix, CopyAssignmentZeroMatrixAddOne)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>> matrix{ builder };

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

TEST(CudaSparseMatrix, CopyAssignmentConstZeroMatrixAddOne)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  const micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>> matrix{ builder };

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

TEST(CudaSparseMatrix, MoveAssignmentConstZeroMatrixAddOne)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>> matrix{ builder };

  auto oneMatrix = std::move(matrix);
  if (matrix.AsDeviceParam().d_data_ != nullptr)
  {
    throw std::runtime_error(
        "The 'd_data_' pointer of oneMatrix is not initialized to a null pointer in the move constructor.");
  }

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

TEST(CudaSparseMatrix, MoveAssignmentZeroMatrixAddOne)
{
  auto builder = micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>>::Create(2)
                     .WithElement(0, 0)
                     .WithElement(1, 1)
                     .InitialValue(0.0);

  micm::CudaSparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseColumn<1>> matrix{ builder };

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
  if (matrix.AsDeviceParam().d_data_ != nullptr)
  {
    throw std::runtime_error(
        "The 'd_data_' pointer of oneMatrix is not initialized to a null pointer in the move constructor.");
  }

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

template <std::size_t cuda_matrix_vector_length>
void TestSingleBlockMatrixAddOneElement()
{
  auto matrix = testSingleBlockMatrix<
      micm::CudaSparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<cuda_matrix_vector_length>>();

  {
    std::size_t elem = matrix.VectorIndex(3, 2);
    EXPECT_EQ(elem, 3 * cuda_matrix_vector_length);
    matrix.AsVector()[elem] = 42;
    EXPECT_EQ(matrix.AsVector()[elem], 42);
  }
  {
    std::size_t elem = matrix.VectorIndex(2, 3);
    EXPECT_EQ(elem, 4 * cuda_matrix_vector_length);
    matrix.AsVector()[elem] = 39;
    EXPECT_EQ(matrix.AsVector()[elem], 39);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), cuda_matrix_vector_length);
  EXPECT_EQ(matrix.GroupSize(), 5 * cuda_matrix_vector_length);
  EXPECT_EQ(matrix.NumberOfGroups(1), 1);

  matrix.CopyToDevice();
  auto param = matrix.AsDeviceParam();
  std::size_t elem_id = 2; // in this example, 2 refers to matrix[2][1] which is non-zero
  micm::cuda::SparseMatrixAddOneElementDriver(param, elem_id, 0, cuda_matrix_vector_length);
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][2][1], 46);
}

TEST(CudaSparseMatrix, SingleBlockMatrixAddOneElement)
{
  TestSingleBlockMatrixAddOneElement<3>();
  TestSingleBlockMatrixAddOneElement<32>();
  TestSingleBlockMatrixAddOneElement<37>();
  TestSingleBlockMatrixAddOneElement<65>();
}

template <std::size_t cuda_matrix_vector_length>
void TestMultiBlockMatrixAddOneElement()
{
  auto matrix = testMultiBlockMatrix<
      micm::CudaSparseMatrix,
      micm::SparseMatrixVectorOrderingCompressedSparseColumn<cuda_matrix_vector_length>>();

  {
    std::size_t elem = matrix.VectorIndex(0, 2, 3);
    EXPECT_EQ(elem, 4 * cuda_matrix_vector_length);
    auto idx = 40;
    elem = matrix.VectorIndex(idx, 2, 3);
    auto a = std::floor(static_cast<double>(idx) / cuda_matrix_vector_length);
    auto b = matrix.VectorIndex(0, 2, 3) / cuda_matrix_vector_length;
    auto c = idx % cuda_matrix_vector_length;
    EXPECT_EQ(elem, 5 * cuda_matrix_vector_length * a + b * cuda_matrix_vector_length + c);
  }
  EXPECT_EQ(matrix.GroupVectorSize(), cuda_matrix_vector_length);
  EXPECT_EQ(matrix.GroupSize(), cuda_matrix_vector_length * 5);
  EXPECT_EQ(matrix.NumberOfGroups(53), std::ceil(53 / static_cast<double>(cuda_matrix_vector_length)));

  matrix.CopyToDevice();
  auto param = matrix.AsDeviceParam();
  std::size_t elem_id = 4; // in this example, 4 refers to matrix[2,3] which is non-zero
  std::size_t grid_id = 33;
  micm::cuda::SparseMatrixAddOneElementDriver(param, elem_id, grid_id, cuda_matrix_vector_length);
  matrix.CopyToHost();

  EXPECT_EQ(matrix[grid_id][2][3], 80);
}

TEST(CudaSparseMatrix, MultiBlockMatrixAddOneElement)
{
  TestMultiBlockMatrixAddOneElement<3>();
  TestMultiBlockMatrixAddOneElement<32>();
  TestMultiBlockMatrixAddOneElement<37>();
  TestMultiBlockMatrixAddOneElement<65>();
}
