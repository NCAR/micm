#include <micm/kokkos/util/kokkos_padded_vector.hpp>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

TEST(KokkosPaddedVector, Compiles) {}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}