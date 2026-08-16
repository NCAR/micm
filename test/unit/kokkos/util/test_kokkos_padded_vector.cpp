#include <micm/kokkos/util/kokkos_padded_vector.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

TEST(KokkosPaddedVector, Compiles)
{
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}