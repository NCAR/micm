#include <gtest/gtest.h>

#include <micm/util/cuda_vector_matrix.cuh>
#include <micm/util/cuda_vector_matrix.hpp>
#include <numeric>

#include "cuda_matrix_utils.cuh"
#include "test_matrix_policy.hpp"

template<class T>
using Group1MatrixAlias = micm::CudaVectorMatrix<T, 1>;
template<class T>
using Group2MatrixAlias = micm::CudaVectorMatrix<T, 2>;
template<class T>
using Group3MatrixAlias = micm::CudaVectorMatrix<T, 3>;
template<class T>
using Group4MatrixAlias = micm::CudaVectorMatrix<T, 4>;

TEST(CudaVectorMatrix, DeviceMemCopy)
{
  std::vector<double> h_vector{ 1, 2, 3, 4 };
  double* h_data = h_vector.data();
  std::size_t num_elements = h_vector.size();
  CudaVectorMatrixParam param;

  micm::cuda::MallocVector(param, num_elements);
  micm::cuda::CopyToDevice(param, h_vector);
  micm::cuda::SquareDriver(param);
  micm::cuda::CopyToHost(param, h_vector);

  EXPECT_EQ(h_vector[0], 1 * 1);
  EXPECT_EQ(h_vector[1], 2 * 2);
  EXPECT_EQ(h_vector[2], 3 * 3);
  EXPECT_EQ(h_vector[3], 4 * 4);
}

TEST(CudaVectorMatrix, IntDataType)
{
  std::vector<std::vector<int>> h_vector{ { 1, 2 }, { 3, 4 } };
  auto matrix = micm::CudaVectorMatrix<int, 2>(h_vector);

  matrix[0][0] = 5;

  EXPECT_EQ(matrix[0][0], 5);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);
}

TEST(CudaVectorMatrix, IntDataTypeCopyAssignment)
{
  std::vector<std::vector<int>> h_vector{ { 1, 2 }, { 3, 4 } };
  auto matrix = micm::CudaVectorMatrix<int, 2>(h_vector);

  matrix[0][0] = 5;

  EXPECT_EQ(matrix[0][0], 5);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  auto matrix2 = matrix;
  matrix2[0][0] = 10;

  EXPECT_EQ(matrix[0][0], 5);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  EXPECT_EQ(matrix2[0][0], 10);
  EXPECT_EQ(matrix2[0][1], 2);
  EXPECT_EQ(matrix2[1][0], 3);
  EXPECT_EQ(matrix2[1][1], 4);
}

TEST(CudaVectorMatrix, IntDataTypeMoveAssignment)
{
  std::vector<std::vector<int>> h_vector{ { 1, 2 }, { 3, 4 } };
  auto matrix = micm::CudaVectorMatrix<int, 2>(h_vector);

  matrix[0][0] = 5;

  EXPECT_EQ(matrix[0][0], 5);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  auto matrix2 = std::move(matrix);
  matrix2[0][0] = 10;

  EXPECT_EQ(matrix2[0][0], 10);
  EXPECT_EQ(matrix2[0][1], 2);
  EXPECT_EQ(matrix2[1][0], 3);
  EXPECT_EQ(matrix2[1][1], 4);
}

template<class T, std::size_t L = MICM_DEFAULT_VECTOR_SIZE>
static void ModifyAndSyncToHost(micm::CudaVectorMatrix<T, L>& matrix)
{
  matrix.CopyToDevice();
  auto matrixParam = matrix.AsDeviceParam();
  micm::cuda::SquareDriver(matrixParam);
  matrix.CopyToHost();
}

TEST(CudaVectorMatrix, CopyConstructorVerifyDeviceMemoryEqual)
{
  std::vector<std::vector<double>> h_vector{ { 1, 2 }, { 3, 4 } };
  auto matrix = micm::CudaVectorMatrix<double, 2>(h_vector);

  matrix[0][0] = 5;

  EXPECT_EQ(matrix[0][0], 5);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  matrix.CopyToDevice();
  auto matrixParam = matrix.AsDeviceParam();
  micm::cuda::SquareDriver(matrixParam);

  EXPECT_EQ(matrix[0][0], 5);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  auto matrix2 = matrix;

  EXPECT_EQ(matrix[0][0], 5);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  EXPECT_EQ(matrix2[0][0], 5);
  EXPECT_EQ(matrix2[0][1], 2);
  EXPECT_EQ(matrix2[1][0], 3);
  EXPECT_EQ(matrix2[1][1], 4);

  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][0], 25);
  EXPECT_EQ(matrix[0][1], 4);
  EXPECT_EQ(matrix[1][0], 9);
  EXPECT_EQ(matrix[1][1], 16);

  EXPECT_EQ(matrix2[0][0], 5);
  EXPECT_EQ(matrix2[0][1], 2);
  EXPECT_EQ(matrix2[1][0], 3);
  EXPECT_EQ(matrix2[1][1], 4);

  matrix2.CopyToHost();

  EXPECT_EQ(matrix[0][0], 25);
  EXPECT_EQ(matrix[0][1], 4);
  EXPECT_EQ(matrix[1][0], 9);
  EXPECT_EQ(matrix[1][1], 16);

  EXPECT_EQ(matrix2[0][0], 25);
  EXPECT_EQ(matrix2[0][1], 4);
  EXPECT_EQ(matrix2[1][0], 9);
  EXPECT_EQ(matrix2[1][1], 16);
}

TEST(CudaVectorMatrix, CopyConstructorSquareAfterCopyAssignment)
{
  std::vector<std::vector<double>> h_vector{ { 1, 2 }, { 3, 4 } };
  auto matrix = micm::CudaVectorMatrix<double, 2>(h_vector);

  auto matrix2 = matrix;

  matrix[0][0] = 5;

  EXPECT_EQ(matrix[0][0], 5);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  EXPECT_EQ(matrix2[0][0], 1);
  EXPECT_EQ(matrix2[0][1], 2);
  EXPECT_EQ(matrix2[1][0], 3);
  EXPECT_EQ(matrix2[1][1], 4);

  ModifyAndSyncToHost(matrix);
  ModifyAndSyncToHost(matrix2);

  EXPECT_EQ(matrix[0][0], 25);
  EXPECT_EQ(matrix[0][1], 4);
  EXPECT_EQ(matrix[1][0], 9);
  EXPECT_EQ(matrix[1][1], 16);
  EXPECT_EQ(matrix2[0][0], 1);
  EXPECT_EQ(matrix2[0][1], 4);
  EXPECT_EQ(matrix2[1][0], 9);
  EXPECT_EQ(matrix2[1][1], 16);
}

TEST(CudaVectorMatrix, CopyConstructorDeSyncedHostDevice)
{
  std::vector<std::vector<double>> h_vector{ { 1, 2 }, { 3, 4 } };
  auto matrix = micm::CudaVectorMatrix<double, 2>(h_vector);

  EXPECT_EQ(matrix[0][0], 1);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  matrix.CopyToDevice();
  auto matrixParam = matrix.AsDeviceParam();
  micm::cuda::SquareDriver(matrixParam);

  EXPECT_EQ(matrix[0][0], 1);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  auto matrix2 = matrix;
  matrix2[0][0] = 5;

  EXPECT_EQ(matrix[0][0], 1);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  EXPECT_EQ(matrix2[0][0], 5);
  EXPECT_EQ(matrix2[0][1], 2);
  EXPECT_EQ(matrix2[1][0], 3);
  EXPECT_EQ(matrix2[1][1], 4);

  matrix.CopyToHost();
  ModifyAndSyncToHost(matrix2);

  EXPECT_EQ(matrix[0][0], 1);
  EXPECT_EQ(matrix[0][1], 4);
  EXPECT_EQ(matrix[1][0], 9);
  EXPECT_EQ(matrix[1][1], 16);

  EXPECT_EQ(matrix2[0][0], 25);
  EXPECT_EQ(matrix2[0][1], 4);
  EXPECT_EQ(matrix2[1][0], 9);
  EXPECT_EQ(matrix2[1][1], 16);
}

TEST(CudaVectorMatrix, CopyAssignment)
{
  std::vector<std::vector<double>> h_vector{ { 1, 2 }, { 3, 4 } };
  auto matrix = micm::CudaVectorMatrix<double, 2>(h_vector);

  micm::CudaVectorMatrix<double, 2> matrix2;
  matrix2 = matrix;

  matrix[0][0] = 5;

  EXPECT_EQ(matrix[0][0], 5);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);
  EXPECT_EQ(matrix2[0][0], 1);
  EXPECT_EQ(matrix2[0][1], 2);
  EXPECT_EQ(matrix2[1][0], 3);
  EXPECT_EQ(matrix2[1][1], 4);

  ModifyAndSyncToHost(matrix);
  ModifyAndSyncToHost(matrix2);

  EXPECT_EQ(matrix[0][0], 25);
  EXPECT_EQ(matrix[0][1], 4);
  EXPECT_EQ(matrix[1][0], 9);
  EXPECT_EQ(matrix[1][1], 16);
  EXPECT_EQ(matrix2[0][0], 1);
  EXPECT_EQ(matrix2[0][1], 4);
  EXPECT_EQ(matrix2[1][0], 9);
  EXPECT_EQ(matrix2[1][1], 16);
}

TEST(CudaVectorMatrix, MoveConstructor)
{
  std::vector<std::vector<double>> h_vector{ { 1, 2 }, { 3, 4 } };
  auto matrix = micm::CudaVectorMatrix<double, 2>(h_vector);

  EXPECT_EQ(matrix[0][0], 1);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  matrix.CopyToDevice();
  auto matrixParam = matrix.AsDeviceParam();
  micm::cuda::SquareDriver(matrixParam);
  matrix[0][0] = 5;

  EXPECT_EQ(matrix[0][0], 5);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[1][0], 3);
  EXPECT_EQ(matrix[1][1], 4);

  auto matrix2 = std::move(matrix);

  EXPECT_EQ(matrix2[0][0], 5);
  EXPECT_EQ(matrix2[0][1], 2);
  EXPECT_EQ(matrix2[1][0], 3);
  EXPECT_EQ(matrix2[1][1], 4);

  matrix2.CopyToHost();

  EXPECT_EQ(matrix2[0][0], 1);
  EXPECT_EQ(matrix2[0][1], 4);
  EXPECT_EQ(matrix2[1][0], 9);
  EXPECT_EQ(matrix2[1][1], 16);
}

TEST(CudaVectorMatrix, MoveAssignment)
{
  std::vector<std::vector<double>> h_vector{ { 1, 2 }, { 3, 4 } };
  auto matrix = micm::CudaVectorMatrix<double, 2>(h_vector);

  micm::CudaVectorMatrix<double, 2> matrix2;
  matrix2 = std::move(matrix);

  EXPECT_EQ(matrix2[0][0], 1);
  EXPECT_EQ(matrix2[0][1], 2);
  EXPECT_EQ(matrix2[1][0], 3);
  EXPECT_EQ(matrix2[1][1], 4);

  ModifyAndSyncToHost(matrix2);

  EXPECT_EQ(matrix2[0][0], 1);
  EXPECT_EQ(matrix2[0][1], 4);
  EXPECT_EQ(matrix2[1][0], 9);
  EXPECT_EQ(matrix2[1][1], 16);
}

TEST(VectorMatrix, SmallVectorMatrix)
{
  auto matrix = testSmallMatrix<Group2MatrixAlias>();

  matrix.CopyToDevice();
  auto devParam = matrix.AsDeviceParam();
  micm::cuda::SquareDriver(devParam);
  matrix.CopyToHost();

  EXPECT_EQ(matrix[1][3], 64.7 * 64.7);
  EXPECT_EQ(matrix[0][0], 41.2 * 41.2);
  EXPECT_EQ(matrix[2][4], 102.3 * 102.3);

  std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 4 * 5);
  EXPECT_EQ(matrix.GroupSize(), 2 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(), 2);
  EXPECT_EQ(matrix.GroupVectorSize(), 2);
  EXPECT_EQ(data[0], 41.2 * 41.2);
  EXPECT_EQ(data[2 * 5 + 0 + 2 * 4], 102.3 * 102.3);
  EXPECT_EQ(data[1 + 2 * 3], 64.7 * 64.7);
}

TEST(CudaVectorMatrix, SmallConstVectorMatrix)
{
  auto matrix = testSmallConstMatrix<Group4MatrixAlias>();

  matrix.CopyToDevice();
  matrix.CopyToHost();

  EXPECT_EQ(matrix[1][3], 64.7);
  EXPECT_EQ(matrix[0][0], 41.2);
  EXPECT_EQ(matrix[2][4], 102.3);

  const std::vector<double>& data = matrix.AsVector();

  EXPECT_EQ(data.size(), 4 * 5);
  EXPECT_EQ(matrix.GroupSize(), 4 * 5);
  EXPECT_EQ(matrix.NumberOfGroups(), 1);
  EXPECT_EQ(matrix.GroupVectorSize(), 4);
  EXPECT_EQ(data[0], 41.2);
  EXPECT_EQ(data[2 + 4 * 4], 102.3);
  EXPECT_EQ(data[1 + 4 * 3], 64.7);
}

TEST(CudaVectorMatrix, InitializeVectorMatrix)
{
  auto matrix = testInializeMatrix<Group1MatrixAlias>();
  matrix.CopyToDevice();
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][0], 12.4);
  EXPECT_EQ(matrix[1][0], 12.4);
  EXPECT_EQ(matrix[1][2], 12.4);
}

TEST(CudaVectorMatrix, InitializeConstVectorMatrix)
{
  auto matrix = testInializeConstMatrix<Group2MatrixAlias>();
  matrix.CopyToDevice();
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][0], 12.4);
  EXPECT_EQ(matrix[1][0], 12.4);
  EXPECT_EQ(matrix[1][2], 12.4);
}

TEST(CudaVectorMatrix, LoopOverVectorMatrix)
{
  Group2MatrixAlias<double> matrix(3, 4, 0);
  for (std::size_t i{}; i < matrix.size(); ++i)
  {
    for (std::size_t j{}; j < matrix[i].size(); ++j)
    {
      matrix[i][j] = i * 100 + j;
    }
  }

  EXPECT_EQ(matrix[0][0], 0);
  EXPECT_EQ(matrix[1][2], 102);
  EXPECT_EQ(matrix[2][3], 203);
  EXPECT_EQ(matrix[0][3], 3);

  matrix.CopyToDevice();
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][0], 0);
  EXPECT_EQ(matrix[1][2], 102);
  EXPECT_EQ(matrix[2][3], 203);
  EXPECT_EQ(matrix[0][3], 3);
}

TEST(CudaVectorMatrix, LoopOverConstVectorMatrix)
{
  Group2MatrixAlias<double> matrix(3, 4, 0);
  for (std::size_t i{}; i < matrix.size(); ++i)
  {
    for (std::size_t j{}; j < matrix[i].size(); ++j)
    {
      matrix[i][j] = i * 100 + j;
    }
  }

  const Group2MatrixAlias<double> const_matrix = matrix;

  EXPECT_EQ(const_matrix[0][0], 0);
  EXPECT_EQ(const_matrix[1][2], 102);
  EXPECT_EQ(const_matrix[2][3], 203);
  EXPECT_EQ(const_matrix[0][3], 3);

  matrix.CopyToDevice();
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][0], 0);
  EXPECT_EQ(matrix[1][2], 102);
  EXPECT_EQ(matrix[2][3], 203);
  EXPECT_EQ(matrix[0][3], 3);
}

TEST(CudaVectorMatrix, ConversionToVector)
{
  auto matrix = testConversionToVector<Group3MatrixAlias>();
  matrix.CopyToDevice();
  matrix.CopyToHost();

  auto slice = matrix[1];

  EXPECT_EQ(slice[0], 13.2);
  EXPECT_EQ(slice[1], 31.2);
  EXPECT_EQ(slice[2], 314.2);
}

TEST(CudaVectorMatrix, ConstConversionToVector)
{
  auto matrix = testConstConversionToVector<Group1MatrixAlias>();
  matrix.CopyToDevice();
  matrix.CopyToHost();

  auto slice = matrix[1];

  EXPECT_EQ(slice[0], 13.2);
  EXPECT_EQ(slice[1], 31.2);
  EXPECT_EQ(slice[2], 314.2);
}

TEST(CudaVectorMatrix, ConversionFromVector)
{
  Group2MatrixAlias<double> zero_matrix = std::vector<std::vector<double>>{};

  EXPECT_EQ(zero_matrix.size(), 0);

  std::vector<std::vector<double>> vec = { { 412.3, 32.4, 41.3 }, { 5.33, -0.3, 31.2 } };

  Group2MatrixAlias<double> matrix = vec;

  EXPECT_EQ(matrix.size(), 2);
  EXPECT_EQ(matrix[0].size(), 3);
  EXPECT_EQ(matrix[0][0], 412.3);
  EXPECT_EQ(matrix[0][1], 32.4);
  EXPECT_EQ(matrix[0][2], 41.3);
  EXPECT_EQ(matrix[1].size(), 3);
  EXPECT_EQ(matrix[1][0], 5.33);
  EXPECT_EQ(matrix[1][1], -0.3);
  EXPECT_EQ(matrix[1][2], 31.2);
}

TEST(CudaVectorMatrix, AssignmentFromVector)
{
  auto matrix = testAssignmentFromVector<Group2MatrixAlias>();
  matrix.CopyToDevice();
  matrix.CopyToHost();

  EXPECT_EQ(matrix[0][0], 0.0);
  EXPECT_EQ(matrix[2][0], 14.3);
  EXPECT_EQ(matrix[2][1], 52.3);
  EXPECT_EQ(matrix[2][2], 65.7);
  EXPECT_EQ(matrix[3][0], 0.0);
}

TEST(CudaVectorMatrix, Axpy)
{
  const double alpha = 2.0;

  // Generate a 20 x 10 matrix with all elements set to 10.0
  auto gpu_x = micm::CudaVectorMatrix<double, 10>(20, 10, 10.0);
  auto gpu_y = micm::CudaVectorMatrix<double, 10>(20, 10, 20.0);
  gpu_x[0][1] = 20.0;
  gpu_x[1][1] = 30.0;

  gpu_x.CopyToDevice();
  gpu_y.CopyToDevice();
  gpu_y.Axpy(alpha, gpu_x, 1, 1);
  gpu_y.CopyToHost();

  EXPECT_EQ(gpu_y[0][0], 40.0);
  EXPECT_EQ(gpu_y[0][1], 60.0);
  EXPECT_EQ(gpu_y[1][1], 80.0);

  double sum = std::accumulate(gpu_y.AsVector().begin(), gpu_y.AsVector().end(), 0);
  EXPECT_EQ(sum, (10 * 2.0 + 20.0) * 198 + 20.0 * 2.0 + 20.0 + 30.0 * 2.0 + 20.0);
}
