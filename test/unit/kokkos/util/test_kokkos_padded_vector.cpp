// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/kokkos/util/kokkos_padded_vector.hpp>
#include <micm/util/types.hpp>
#include <micm/util/view_category.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <type_traits>
#include <utility>
#include <vector>

// KokkosPaddedVector is not a wrapper around std::vector: it holds a host std::vector, an
// unmanaged host View aliasing that vector, and a *separately allocated* Kokkos::View for the
// device. Those two allocations stay distinct even in a Serial/host-only build, so the
// round-trip tests below are meaningful on a CPU-only CI runner and not just on a GPU.
//
// This type is also not incidental to the backend: State::conditions_ is a
// KokkosPaddedVector<Conditions, L>.

// Anything that launches a kernel lives in a free function rather than directly in the TEST
// body. gtest expands TEST() into a class whose TestBody() is a *private* member function, and
// nvcc rejects an extended __host__ __device__ lambda -- which is what KOKKOS_LAMBDA becomes
// under CUDA -- inside a private or protected member function. A host-only Kokkos build
// compiles either form, so this only shows up on the CUDA runner. test_kokkos_dense_matrix.cpp
// sidesteps the same restriction by delegating to the shared test_matrix_policy.hpp helpers.

using IntVec4 = micm::KokkosPaddedVector<int, 4>;

static_assert(micm::PaddedVectorLike<IntVec4>, "KokkosPaddedVector must satisfy PaddedVectorLike");
static_assert(
    std::is_same_v<typename IntVec4::category, micm::PaddedVectorTag>,
    "KokkosPaddedVector must be tagged as a padded vector");
static_assert(std::is_same_v<typename IntVec4::value_type, int>, "value_type must be the element type");
static_assert(IntVec4::GroupVectorSize() == 4, "GroupVectorSize must report the template vector length");

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

TEST(KokkosPaddedVector, DefaultConstructor)
{
  IntVec4 vec;
  EXPECT_EQ(vec.size(), 0);
  EXPECT_EQ(vec.PaddedSize(), 0);
}

TEST(KokkosPaddedVector, PaddedSizeRounding)
{
  // PaddedSize is ceil(n / L) * L, including the exact-multiple and empty cases.
  EXPECT_EQ((micm::KokkosPaddedVector<int, 1>(7).PaddedSize()), 7);
  EXPECT_EQ((micm::KokkosPaddedVector<int, 3>(7).PaddedSize()), 9);
  EXPECT_EQ((micm::KokkosPaddedVector<int, 5>(10).PaddedSize()), 10);
  EXPECT_EQ((micm::KokkosPaddedVector<int, 8>(1).PaddedSize()), 8);
  EXPECT_EQ((micm::KokkosPaddedVector<int, 8>(0).PaddedSize()), 0);
  // size() always reports the logical length, never the padded one.
  EXPECT_EQ((micm::KokkosPaddedVector<int, 8>(1).size()), 1);
}

TEST(KokkosPaddedVector, SizedConstructorFillsPadding)
{
  IntVec4 vec(3, 7);
  EXPECT_EQ(vec.size(), 3);
  EXPECT_EQ(vec.PaddedSize(), 4);
  // operator[] reaches the padding as well as the logical elements.
  for (micm::Index i = 0; i < 4; ++i)
  {
    EXPECT_EQ(vec[i], 7) << "element " << i;
  }
}

TEST(KokkosPaddedVector, IterationStopsAtLogicalSize)
{
  IntVec4 vec(3, 7);
  vec[3] = 99;  // padding element
  std::vector<int> seen;
  for (const auto elem : vec)
  {
    seen.push_back(elem);
  }
  EXPECT_EQ(seen, (std::vector<int>{ 7, 7, 7 })) << "iteration must not expose the padding";
}

TEST(KokkosPaddedVector, VectorConstructor)
{
  IntVec4 vec(std::vector<int>{ 1, 2, 3, 4, 5 });
  EXPECT_EQ(vec.size(), 5);
  EXPECT_EQ(vec.PaddedSize(), 8);
  EXPECT_EQ(vec[0], 1);
  EXPECT_EQ(vec[4], 5);
  EXPECT_EQ(vec[5], 0) << "padding must be value-initialized";
}

TEST(KokkosPaddedVector, CopyToDeviceAndHost)
{
  IntVec4 vec(3, 0);
  vec[0] = 1;
  vec[1] = 2;
  vec[2] = 3;

  vec.CopyToDevice();

  // Clear the host buffer so a no-op CopyToHost cannot be mistaken for a successful download.
  vec[0] = 0;
  vec[1] = 0;
  vec[2] = 0;

  vec.CopyToHost();
  EXPECT_EQ(vec[0], 1);
  EXPECT_EQ(vec[1], 2);
  EXPECT_EQ(vec[2], 3);
}

void CheckGetViewAddressesDeviceMemory()
{
  IntVec4 vec(3, 0);
  vec[0] = 1;
  vec[1] = 2;
  vec[2] = 3;
  vec.CopyToDevice();

  // Double every element through the device view.
  auto view = vec.GetView();
  Kokkos::parallel_for(
      "double_elements",
      3,
      KOKKOS_LAMBDA(const micm::Index i) { view[i] = view[i] * 2; });
  Kokkos::fence();

  // The host buffer is a separate allocation, so it must still hold the original values.
  EXPECT_EQ(vec[0], 1) << "GetView() must not alias the host buffer";
  EXPECT_EQ(vec[1], 2);
  EXPECT_EQ(vec[2], 3);

  vec.CopyToHost();
  EXPECT_EQ(vec[0], 2);
  EXPECT_EQ(vec[1], 4);
  EXPECT_EQ(vec[2], 6);
}

TEST(KokkosPaddedVector, GetViewAddressesDeviceMemory)
{
  CheckGetViewAddressesDeviceMemory();
}

void CheckDeviceViewReportsLogicalSize()
{
  IntVec4 vec(3, 0);
  vec.CopyToDevice();

  // size() must be the logical length even though the underlying view spans the padding.
  auto view = vec.GetView();
  EXPECT_EQ(view.size(), 3);
  EXPECT_EQ(vec.PaddedSize(), 4);

  // ...and it must report the same value inside a kernel.
  Kokkos::View<micm::Index*> observed("observed", 1);
  Kokkos::parallel_for(
      "read_size",
      1,
      KOKKOS_LAMBDA(const micm::Index) { observed(0) = view.size(); });
  Kokkos::fence();
  auto observed_host = Kokkos::create_mirror_view(observed);
  Kokkos::deep_copy(observed_host, observed);
  EXPECT_EQ(observed_host(0), 3);
}

TEST(KokkosPaddedVector, DeviceViewReportsLogicalSize)
{
  CheckDeviceViewReportsLogicalSize();
}

void CheckConstDeviceViewConversionKeepsSize()
{
  IntVec4 vec(3, 0);
  vec[0] = 5;
  vec.CopyToDevice();

  // State::VectorView is the ConstViewType, so this conversion is on the path any solver
  // takes when it hands a padded vector to a kernel as read-only data. It has to carry the
  // logical size across, not just the underlying view.
  typename IntVec4::ViewType mutable_view = vec.GetView();
  typename IntVec4::ConstViewType const_view = mutable_view;

  // size(), begin() and end() only read the size member and the view's pointer, so they are
  // safe to call on the host. Dereferencing the view is not: under CUDA it addresses
  // CudaSpace, so the element has to be read inside a kernel.
  EXPECT_EQ(const_view.size(), mutable_view.size());
  EXPECT_EQ(const_view.size(), 3);
  EXPECT_EQ(const_view.end() - const_view.begin(), 3);

  Kokkos::View<int*> observed("observed", 1);
  Kokkos::parallel_for(
      "read_const_view",
      1,
      KOKKOS_LAMBDA(const micm::Index) { observed(0) = const_view[0]; });
  Kokkos::fence();
  auto observed_host = Kokkos::create_mirror_view(observed);
  Kokkos::deep_copy(observed_host, observed);
  EXPECT_EQ(observed_host(0), 5);
}

TEST(KokkosPaddedVector, ConstDeviceViewConversionKeepsSize)
{
  CheckConstDeviceViewConversionKeepsSize();
}

void CheckCopyConstructorCopiesDeviceData()
{
  IntVec4 vec(3, 0);
  vec.CopyToDevice();

  // Put a value on the device that does not exist in the host buffer, so the assertion below
  // can only pass if the copy constructor deep-copied the device view.
  auto view = vec.GetView();
  Kokkos::parallel_for(
      "seed_device",
      3,
      KOKKOS_LAMBDA(const micm::Index i) { view[i] = 40 + static_cast<int>(i); });
  Kokkos::fence();

  IntVec4 copy(vec);
  EXPECT_EQ(copy.size(), 3);
  EXPECT_EQ(copy[0], 0) << "the copied host buffer is the stale one, as in the source";

  copy.CopyToHost();
  EXPECT_EQ(copy[0], 40);
  EXPECT_EQ(copy[1], 41);
  EXPECT_EQ(copy[2], 42);
}

TEST(KokkosPaddedVector, CopyConstructorCopiesDeviceData)
{
  CheckCopyConstructorCopiesDeviceData();
}

void CheckCopyAssignmentCopiesDeviceData()
{
  IntVec4 vec(3, 0);
  vec.CopyToDevice();
  auto view = vec.GetView();
  Kokkos::parallel_for(
      "seed_device",
      3,
      KOKKOS_LAMBDA(const micm::Index i) { view[i] = 70 + static_cast<int>(i); });
  Kokkos::fence();

  IntVec4 other(1, 0);
  other = vec;
  ASSERT_EQ(other.size(), 3) << "assignment must resize";
  other.CopyToHost();
  EXPECT_EQ(other[0], 70);
  EXPECT_EQ(other[1], 71);
  EXPECT_EQ(other[2], 72);
}

TEST(KokkosPaddedVector, CopyAssignmentCopiesDeviceData)
{
  CheckCopyAssignmentCopiesDeviceData();
}

TEST(KokkosPaddedVector, MoveConstructorPreservesValues)
{
  IntVec4 source(3, 0);
  source[0] = 11;
  source[1] = 12;
  source[2] = 13;

  IntVec4 moved(std::move(source));
  EXPECT_EQ(moved.size(), 3);
  EXPECT_EQ(moved[0], 11);
  EXPECT_EQ(moved[1], 12);
  EXPECT_EQ(moved[2], 13);

  // The unmanaged host view must still address the buffer the vector now owns.
  moved.CopyToDevice();
  moved[0] = 0;
  moved.CopyToHost();
  EXPECT_EQ(moved[0], 11);
}

TEST(KokkosPaddedVector, Equality)
{
  IntVec4 a(3, 0);
  a[0] = 1;
  a[1] = 2;
  a[2] = 3;
  IntVec4 b(a);
  EXPECT_TRUE(a == b);

  b[1] = 9;
  EXPECT_FALSE(a == b);

  // Comparison against std::vector is over the padded buffer, so include the padding.
  const std::vector<int> expected{ 1, 2, 3, 0 };
  EXPECT_TRUE(a == expected);
  EXPECT_TRUE(expected == a);
  EXPECT_FALSE(a == (std::vector<int>{ 1, 2, 3 })) << "a shorter vector must not compare equal";
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
