#include <gtest/gtest.h>

#include <micm/process/arrhenius_rate_constant.hpp>
#include <micm/solver/cuda_rosenbrock.cuh>
#include <micm/solver/cuda_rosenbrock.hpp>
#include <micm/solver/rosenbrock.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>
#include <micm/util/cuda_vector_matrix.hpp>

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

// the following alias works for a CudaVectorMatrix with 1 row and any columns
template<class T>
using Group1CudaVectorMatrix = micm::CudaVectorMatrix<T, 1>;
// the following alias works for a CudaVectorMatrix with 7 rows and any columns
template<class T>
using Group7CudaVectorMatrix = micm::CudaVectorMatrix<T, 7>;

template<
    template<class>
    class MatrixPolicy,
    template<class>
    class SparseMatrixPolicy,
    class LinearSolverPolicy,
    class RosenbrockPolicy>
RosenbrockPolicy getSolver(std::size_t number_of_grid_cells)
{
  // ---- foo  bar  baz  quz  quuz
  // foo   0    1    2    -    -
  // bar   3    4    5    -    -
  // baz   6    -    7    -    -
  // quz   -    8    -    9    -
  // quuz 10    -   11    -    12

  auto foo = micm::Species("foo");
  auto bar = micm::Species("bar");
  auto baz = micm::Species("baz");
  auto quz = micm::Species("quz");
  auto quuz = micm::Species("quuz");

  micm::Phase gas_phase{ std::vector<micm::Species>{ foo, bar, baz, quz, quuz } };

  micm::Process r1 = micm::Process::create()
                         .reactants({ foo, baz })
                         .products({ yields(bar, 1), yields(quuz, 2.4) })
                         .phase(gas_phase)
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 2.0e-11, .B_ = 0, .C_ = 110 }));

  micm::Process r2 = micm::Process::create()
                         .reactants({ bar })
                         .products({ yields(foo, 1), yields(quz, 1.4) })
                         .phase(gas_phase)
                         .rate_constant(micm::ArrheniusRateConstant({ .A_ = 1.0e-6 }));

  micm::Process r3 = micm::Process::create().reactants({ quz }).products({}).phase(gas_phase).rate_constant(
      micm::ArrheniusRateConstant({ .A_ = 3.5e-6 }));

  return RosenbrockPolicy(
      micm::System(micm::SystemParameters{ .gas_phase_ = gas_phase }),
      std::vector<micm::Process>{ r1, r2, r3 },
      micm::RosenbrockSolverParameters::three_stage_rosenbrock_parameters(number_of_grid_cells, false));
}

template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
void testAlphaMinusJacobian(std::size_t number_of_grid_cells)
{
  auto gpu_solver = getSolver<
      MatrixPolicy,
      SparseMatrixPolicy,
      LinearSolverPolicy,
      micm::CudaRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>>(number_of_grid_cells);
  auto jacobian = gpu_solver.GetState().jacobian_;
  EXPECT_EQ(jacobian.size(), number_of_grid_cells);
  EXPECT_EQ(jacobian[0].size(), 5);
  EXPECT_EQ(jacobian[0][0].size(), 5);
  EXPECT_GE(jacobian.AsVector().size(), 13 * number_of_grid_cells);
  for (auto& elem : jacobian.AsVector())
    elem = 100.0;
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    jacobian[i_cell][0][0] = 12.2;
    jacobian[i_cell][0][1] = 24.3 * (i_cell + 2);
    jacobian[i_cell][0][2] = 42.3;
    jacobian[i_cell][1][0] = 0.43;
    jacobian[i_cell][1][1] = 23.4;
    jacobian[i_cell][1][2] = 83.4 / (i_cell + 3);
    jacobian[i_cell][2][0] = 4.74;
    jacobian[i_cell][2][2] = 6.91;
    jacobian[i_cell][3][1] = 59.1;
    jacobian[i_cell][3][3] = 83.4;
    jacobian[i_cell][4][0] = 78.5;
    jacobian[i_cell][4][2] = 53.6;
    jacobian[i_cell][4][4] = 1.0;
  }
  auto cpu_jacobian = jacobian;

  gpu_solver.AlphaMinusJacobian(jacobian, 42.042);
  for (std::size_t i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    EXPECT_EQ(jacobian[i_cell][0][0], 42.042 - 12.2);
    EXPECT_EQ(jacobian[i_cell][0][1], -24.3 * (i_cell + 2));
    EXPECT_EQ(jacobian[i_cell][0][2], -42.3);
    EXPECT_EQ(jacobian[i_cell][1][0], -0.43);
    EXPECT_EQ(jacobian[i_cell][1][1], 42.042 - 23.4);
    EXPECT_EQ(jacobian[i_cell][1][2], -83.4 / (i_cell + 3));
    EXPECT_EQ(jacobian[i_cell][2][0], -4.74);
    EXPECT_EQ(jacobian[i_cell][2][2], 42.042 - 6.91);
    EXPECT_EQ(jacobian[i_cell][3][1], -59.1);
    EXPECT_EQ(jacobian[i_cell][3][3], 42.042 - 83.4);
    EXPECT_EQ(jacobian[i_cell][4][0], -78.5);
    EXPECT_EQ(jacobian[i_cell][4][2], -53.6);
    EXPECT_EQ(jacobian[i_cell][4][4], 42.042 - 1.0);
  }

  auto cpu_solver = getSolver<
      MatrixPolicy,
      SparseMatrixPolicy,
      micm::LinearSolver<double, SparseMatrixPolicy>,
      micm::RosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, micm::LinearSolver<double, SparseMatrixPolicy>>>(
      number_of_grid_cells);
  cpu_solver.AlphaMinusJacobian(cpu_jacobian, 42.042);
  std::vector<double> jacobian_gpu_vector = jacobian.AsVector();
  std::vector<double> jacobian_cpu_vector = cpu_jacobian.AsVector();
  for (int i = 0; i < jacobian_cpu_vector.size(); i++)
  {
    EXPECT_EQ(jacobian_cpu_vector[i], jacobian_gpu_vector[i]);
  }
}

// In this test, all the elements in the same array are identical;
// thus the calculated RMSE should be the same no matter what the size of the array is. 
template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
double testNormalizedErrorConst(const size_t nrows, const size_t ncols)
{
  auto gpu_solver = getSolver<
      MatrixPolicy,
      SparseMatrixPolicy,
      LinearSolverPolicy,
      micm::CudaRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>>(nrows * ncols);
  
  double atol = gpu_solver.parameters_.absolute_tolerance_;
  double rtol = gpu_solver.parameters_.relative_tolerance_;

  auto y_old  = MatrixPolicy<double>(nrows,ncols,1.0);
  auto y_new  = MatrixPolicy<double>(nrows,ncols,2.0);
  auto errors = MatrixPolicy<double>(nrows,ncols,3.0);

  y_old.CopyToDevice();
  y_new.CopyToDevice();
  errors.CopyToDevice();

  double error = gpu_solver.NormalizedError(y_old, y_new, errors);

  double denom = atol+rtol*2.0;
  // use the following function instead to avoid tiny numerical differece
  EXPECT_DOUBLE_EQ( error, std::sqrt(3.0*3.0/(denom*denom)) );
}

// In this test, the elements in the same array are different;
// thus the calculated RMSE will change when the size of the array changes. 
template<template<class> class MatrixPolicy, template<class> class SparseMatrixPolicy, class LinearSolverPolicy>
double testNormalizedErrorDiff(const size_t nrows, const size_t ncols)
{
  auto gpu_solver = getSolver<
      MatrixPolicy,
      SparseMatrixPolicy,
      LinearSolverPolicy,
      micm::CudaRosenbrockSolver<MatrixPolicy, SparseMatrixPolicy, LinearSolverPolicy>>(nrows * ncols);
  
  double atol = gpu_solver.parameters_.absolute_tolerance_;
  double rtol = gpu_solver.parameters_.relative_tolerance_;

  auto y_old  = MatrixPolicy<double>(nrows,ncols,7.7);
  auto y_new  = MatrixPolicy<double>(nrows,ncols,-13.9);
  auto errors = MatrixPolicy<double>(nrows,ncols,81.57);

  double expected_error = 0.0;
  for (size_t i = 0; i < nrows; ++i)
  {
    for (size_t j = 0; j < ncols; ++j)
    {
      y_old[i,j] = y_old[i,j] * i + j;
      y_new[i,j] = y_old[i,j] / (j+1) - i;
      errors[i,j] = errors[i,j] / (i+7) / (j+3);
      double ymax = std::max(std::abs(y_old[i,j]), std::abs(y_new[i,j]));
      double scale = atol + rtol * ymax;
      expected_error += std::pow(errors[i,j] / scale, 2);
    }
  }

  y_old.CopyToDevice();
  y_new.CopyToDevice();
  errors.CopyToDevice();

  double computed_error = gpu_solver.NormalizedError(y_old, y_new, errors);

  EXPECT_DOUBLE_EQ( computed_error, expected_error );
}

TEST(RosenbrockSolver, DenseAlphaMinusJacobian)
{
  testAlphaMinusJacobian<
      Group1VectorMatrix,
      Group1SparseVectorMatrix,
      micm::CudaLinearSolver<double, Group1SparseVectorMatrix>>(1);
  testAlphaMinusJacobian<
      Group2VectorMatrix,
      Group2SparseVectorMatrix,
      micm::CudaLinearSolver<double, Group2SparseVectorMatrix>>(2);
  testAlphaMinusJacobian<
      Group3VectorMatrix,
      Group3SparseVectorMatrix,
      micm::CudaLinearSolver<double, Group3SparseVectorMatrix>>(3);
  testAlphaMinusJacobian<
      Group4VectorMatrix,
      Group4SparseVectorMatrix,
      micm::CudaLinearSolver<double, Group4SparseVectorMatrix>>(4);
}

TEST(RosenbrockSolver, CudaNormalizedError)
{
  std::vector<int> row_array = {1, 7};
  // Based on my experience, assuming we use BLOCK_SIZE = N, it is better to test the length size (L)
  // of arrays at least within the following ranges for robustness: [L<N, N<L<2N, 2N<L<4N, 4N<L<N^2, N^2<L<N^3];
  // Trying some odd and weird numbers is always helpful to reveal a potential bug.
  std::vector<int> col_array = {1, 7, 23, 35, 63, 79, 101, 27997, 33017, 1000201, 2109377, 16975213};

  // tests where RMSE does not change with the size of the array
  for (auto row : row_array)
  {
    for (auto col : col_array)
    {
      testNormalizedErrorConst<
          Group1CudaVectorMatrix,
          Group1SparseVectorMatrix,
          micm::CudaLinearSolver<double, Group1SparseVectorMatrix>>(row,col);
    }
  }
  // tests where RMSE changes with the size of the array
  for (auto row : row_array)
  {
    for (auto col : col_array)
    {
      testNormalizedErrorConst<
          Group1CudaVectorMatrix,
          Group1SparseVectorMatrix,
          micm::CudaLinearSolver<double, Group1SparseVectorMatrix>>(row,col);
    }
  }  
}
