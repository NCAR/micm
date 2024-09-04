#include <micm/cuda/util/cuda_util.cuh>

#include <gtest/gtest.h>

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  micm::cuda::CudaStreamSingleton::GetInstance().CleanUp();  // Ensure cleanup
  return result;
}