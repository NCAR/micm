// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <micm/kokkos/util/kokkos_scalar_view.hpp>
#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <type_traits>

// KokkosScalarView backs a single value with two distinct allocations -- a HostSpace View and
// a default-space View -- and dispatches its accessors with KOKKOS_IF_ON_HOST /
// KOKKOS_IF_ON_DEVICE. Because the two allocations are separate even in a Serial build, the
// round-trip tests here exercise a real copy on a CPU-only runner.
//
// It is the ScalarType of both KokkosDenseMatrix and KokkosSparseMatrix, and the target that
// KokkosSum / KokkosMax / KokkosLOr / KokkosLAnd reduce into.

// Anything that launches a kernel lives in a free function rather than directly in the TEST
// body. gtest expands TEST() into a class whose TestBody() is a *private* member function, and
// nvcc rejects an extended __host__ __device__ lambda -- which is what KOKKOS_LAMBDA becomes
// under CUDA -- inside a private or protected member function. A host-only Kokkos build
// compiles either form, so this only shows up on the CUDA runner.

static_assert(
    std::is_same_v<typename micm::KokkosScalarView<micm::Real>::value_type, micm::Real>,
    "value_type must be the element type");

TEST(KokkosScalarView, DefaultConstructorIsValueInitialized)
{
  micm::KokkosScalarView<micm::Real> scalar;
  EXPECT_EQ(scalar.HostValue(), micm::Real{ 0.0 });
  EXPECT_EQ(static_cast<micm::Real>(scalar), micm::Real{ 0.0 });

  micm::KokkosScalarView<int> int_scalar;
  EXPECT_EQ(int_scalar.HostValue(), 0);
}

TEST(KokkosScalarView, InitialisingConstructor)
{
  micm::KokkosScalarView<micm::Real> scalar(3.5);
  EXPECT_EQ(scalar.HostValue(), micm::Real{ 3.5 });
  // operator T() resolves to the host branch when called from host code.
  EXPECT_EQ(static_cast<micm::Real>(scalar), micm::Real{ 3.5 });
}

TEST(KokkosScalarView, HostDataPointer)
{
  micm::KokkosScalarView<micm::Real> scalar(2.25);
  EXPECT_EQ(*scalar.data(), micm::Real{ 2.25 });
  // HostValue() is a reference into the same host allocation data() returns.
  EXPECT_EQ(scalar.data(), &scalar.HostValue());

  scalar.HostValue() = 9.0;
  EXPECT_EQ(*scalar.data(), micm::Real{ 9.0 });
}

TEST(KokkosScalarView, ConstDataPointer)
{
  const micm::KokkosScalarView<micm::Real> scalar(1.5);
  EXPECT_EQ(*scalar.data(), micm::Real{ 1.5 });
}

TEST(KokkosScalarView, GetDeviceViewIsRankOneSizeOne)
{
  micm::KokkosScalarView<micm::Real> scalar(1.0);
  auto device_view = scalar.GetDeviceView();
  EXPECT_EQ(device_view.extent(0), 1u);
  EXPECT_EQ(device_view.size(), 1u);
}

TEST(KokkosScalarView, CopyToDeviceAndHost)
{
  micm::KokkosScalarView<micm::Real> scalar(5.0);
  scalar.CopyToDevice();

  // Clear the host side so a no-op download cannot masquerade as a successful one.
  scalar.HostValue() = 0.0;
  ASSERT_EQ(scalar.HostValue(), micm::Real{ 0.0 });

  scalar.CopyToHost();
  EXPECT_EQ(scalar.HostValue(), micm::Real{ 5.0 });
}

void CheckDeviceWriteIsVisibleAfterCopyToHost()
{
  micm::KokkosScalarView<micm::Real> scalar(0.0);
  scalar.CopyToDevice();

  auto device_view = scalar.GetDeviceView();
  Kokkos::parallel_for(
      "write_scalar",
      1,
      KOKKOS_LAMBDA(const micm::Index) { device_view(0) = micm::Real{ 42.0 }; });
  Kokkos::fence();

  // The host allocation is distinct, so it must not have moved yet.
  EXPECT_EQ(scalar.HostValue(), micm::Real{ 0.0 }) << "the device view must not alias the host value";

  scalar.CopyToHost();
  EXPECT_EQ(scalar.HostValue(), micm::Real{ 42.0 });
}

TEST(KokkosScalarView, DeviceWriteIsVisibleAfterCopyToHost)
{
  CheckDeviceWriteIsVisibleAfterCopyToHost();
}

void CheckDeviceReadSeesUploadedValue()
{
  micm::KokkosScalarView<micm::Real> scalar(7.25);
  scalar.CopyToDevice();

  auto device_view = scalar.GetDeviceView();
  Kokkos::View<micm::Real*> observed("observed", 1);
  Kokkos::parallel_for(
      "read_scalar",
      1,
      KOKKOS_LAMBDA(const micm::Index) { observed(0) = device_view(0); });
  Kokkos::fence();

  auto observed_host = Kokkos::create_mirror_view(observed);
  Kokkos::deep_copy(observed_host, observed);
  EXPECT_EQ(observed_host(0), micm::Real{ 7.25 });
}

TEST(KokkosScalarView, DeviceReadSeesUploadedValue)
{
  CheckDeviceReadSeesUploadedValue();
}

void CheckCopyAssignmentCopiesTheHostValue()
{
  micm::KokkosScalarView<micm::Real> source(11.0);
  micm::KokkosScalarView<micm::Real> target(0.0);

  target = source;
  EXPECT_EQ(target.HostValue(), micm::Real{ 11.0 });

  // Self-assignment is guarded and must leave the value alone.
  target = target;
  EXPECT_EQ(target.HostValue(), micm::Real{ 11.0 });

  // Note: this operator copies only the host allocation. A caller that needs the assigned
  // value on the device has to call CopyToDevice() afterwards -- asserted here so the
  // requirement is documented rather than discovered.
  target.CopyToDevice();
  auto device_view = target.GetDeviceView();
  Kokkos::View<micm::Real*> observed("observed", 1);
  Kokkos::parallel_for(
      "read_after_assign",
      1,
      KOKKOS_LAMBDA(const micm::Index) { observed(0) = device_view(0); });
  Kokkos::fence();
  auto observed_host = Kokkos::create_mirror_view(observed);
  Kokkos::deep_copy(observed_host, observed);
  EXPECT_EQ(observed_host(0), micm::Real{ 11.0 });
}

TEST(KokkosScalarView, CopyAssignmentCopiesTheHostValue)
{
  CheckCopyAssignmentCopiesTheHostValue();
}

TEST(KokkosScalarView, IntegerSpecialisationRoundTrips)
{
  micm::KokkosScalarView<int> scalar(3);
  scalar.CopyToDevice();
  scalar.HostValue() = 0;
  scalar.CopyToHost();
  EXPECT_EQ(scalar.HostValue(), 3);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
