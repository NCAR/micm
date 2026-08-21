#include <micm/kokkos/util/kokkos_padded_vector.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

TEST(KokkosPaddedVector, DimensionsAndValues)
{
  micm::KokkosPaddedVector<int, 5> vec{ 6, 5, 4, 3, 2, 1 };
  ASSERT_EQ(vec.size(), 6);
  ASSERT_EQ(vec.PaddedSize(), 10);
  int comp = 7;
  for (const auto elem : vec)
  {
    ASSERT_EQ(elem, --comp);
  }
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}