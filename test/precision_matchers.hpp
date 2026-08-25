#pragma once

#include <micm/util/types.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

// Precision-aware floating-point comparisons for tests that must run under both
// MICM_USE_DOUBLE=ON (micm::Real = double) and MICM_USE_DOUBLE=OFF (micm::Real = float).
//
// Every macro here widens to a fixed tolerance in the double build, so a test that passes at
// double precision today keeps exactly the tolerance it had; only the float build is relaxed.

namespace micm::test
{
  /// @brief Machine epsilon of the active working precision (2.2e-16 double, 1.2e-7 float).
  inline constexpr Real REAL_EPS = std::numeric_limits<Real>::epsilon();

  /// @brief Relative-error floor a value reconstructed by arithmetic can be held to.
  ///
  /// Zero in double so double-precision tests are never loosened; a few epsilon in float, which is
  /// the accuracy floor for any quantity reconstructed by a handful of arithmetic operations on
  /// floats. This is *not* the right floor for a value carried through a many-step integration --
  /// use REAL_SOLVE_REL_FLOOR for those.
  inline constexpr double REAL_REL_FLOOR = std::is_same_v<Real, double> ? 0.0 : 8.0 * (double)REAL_EPS;

  /// @brief Relative-error floor between two evaluations of the same closed-form expression.
  ///
  /// A test that re-derives a rate constant from its documented formula does not reproduce micm's
  /// own evaluation bit-for-bit: the two associate the operations differently, and an overload like
  /// std::pow(float, int) silently widens part of the test's expression to double. The gap is a
  /// small multiple of epsilon rather than the 4 ULP an equality macro allows -- micm's
  /// CalculateTunneling uses T*T*T where the test uses std::pow(T, 3), which differ by 8 ULP at
  /// T = 254.7 K. Zero in double.
  inline constexpr double REAL_FORMULA_REL_FLOOR = std::is_same_v<Real, double> ? 0.0 : 32.0 * (double)REAL_EPS;

  /// @brief Relative-error floor for a value produced by integrating the system over many steps.
  ///
  /// Bounded below by the solver's own default relative tolerance at the active precision -- 1e-5
  /// in float, see micm::StateParameters -- times room for global error to accumulate over the
  /// ~1e3 steps these tests take. Asking for less than this asks the solver for accuracy it was
  /// never configured to deliver. The worst case observed across the analytical examples is
  /// 8.8e-5 (AnalyticalExamples.TernaryChemicalActivation, species A). Zero in double.
  inline constexpr double REAL_SOLVE_REL_FLOOR = std::is_same_v<Real, double> ? 0.0 : 2.0e-4;

  /// @brief Combined absolute + relative tolerance for comparing a model value to a reference.
  inline double CombinedTolerance(double abs_tol, double reference)
  {
    return abs_tol + REAL_REL_FLOOR * std::abs(reference);
  }

  /// @brief The larger of a caller-supplied relative tolerance and the precision floor.
  inline double RelativeTolerance(double rel_tol)
  {
    return std::max(rel_tol, REAL_REL_FLOOR);
  }

  /// @brief As CombinedTolerance, but floored at the accuracy an integrated result can reach.
  inline double SolveTolerance(double abs_tol, double reference)
  {
    return abs_tol + REAL_SOLVE_REL_FLOOR * std::abs(reference);
  }

  /// @brief As RelativeTolerance, but floored at the accuracy an integrated result can reach.
  inline double SolveRelativeTolerance(double rel_tol)
  {
    return std::max(rel_tol, REAL_SOLVE_REL_FLOOR);
  }

  /// @brief Tolerance for two evaluations of the same closed-form expression.
  inline double FormulaTolerance(double reference)
  {
    return REAL_FORMULA_REL_FLOOR * std::abs(reference);
  }
}  // namespace micm::test

// Equality to within 4 ULP of the active precision.
//   double build -> EXPECT_DOUBLE_EQ / ASSERT_DOUBLE_EQ
//   float  build -> EXPECT_FLOAT_EQ  / ASSERT_FLOAT_EQ
//
// Use in place of EXPECT_DOUBLE_EQ whenever comparing micm::Real values against a reference. Note:
// like the gtest macros it wraps, this is an *equality* check and is not suitable for values
// legitimately near zero or for reconstruction/backward-error checks whose error exceeds a few ULP
// -- use EXPECT_REAL_CLOSE or EXPECT_REAL_REL for those.
#ifdef MICM_USE_SINGLE
  #define EXPECT_REAL_EQ(a, b) EXPECT_FLOAT_EQ((a), (b))
  #define ASSERT_REAL_EQ(a, b) ASSERT_FLOAT_EQ((a), (b))
#else
  #define EXPECT_REAL_EQ(a, b) EXPECT_DOUBLE_EQ((a), (b))
  #define ASSERT_REAL_EQ(a, b) ASSERT_DOUBLE_EQ((a), (b))
#endif

// |actual - reference| <= abs_tol + (precision floor) * |reference|
//
// For comparing a solver result against a reference whose magnitude varies, where a bare absolute
// tolerance would demand accuracy finer than the working precision can represent. In a double build
// the relative term is zero, so this reduces exactly to EXPECT_NEAR(actual, reference, abs_tol).
#define EXPECT_REAL_CLOSE(actual, reference, abs_tol) \
  EXPECT_NEAR((actual), (reference), micm::test::CombinedTolerance((abs_tol), (reference)))

// |actual - reference| <= max(rel_tol, precision floor) * |reference|
//
// For assertions already expressed as a relative error, where the requested rel_tol may itself sit
// below machine epsilon of the working precision.
#define EXPECT_REAL_REL(actual, reference, rel_tol) \
  EXPECT_NEAR((actual), (reference), micm::test::RelativeTolerance(rel_tol) * std::abs((double)(reference)))

// |actual - reference| <= abs_tol + (integration floor) * |reference|
//
// For comparing a concentration the solver integrated over many steps against a closed-form
// reference. Distinct from EXPECT_REAL_CLOSE because the achievable accuracy is set by the solver's
// own relative tolerance and the global error accumulated over the run, not by the roundoff of a
// few arithmetic operations. In a double build the relative term is zero, so this reduces exactly
// to EXPECT_NEAR(actual, reference, abs_tol).
#define EXPECT_REAL_SOLVE_CLOSE(actual, reference, abs_tol) \
  EXPECT_NEAR((actual), (reference), micm::test::SolveTolerance((abs_tol), (reference)))

// |actual - reference| <= max(rel_tol, integration floor) * |reference|
//
// The EXPECT_REAL_REL counterpart for a quantity carried through many solver steps -- a total the
// kinetics conserve, say -- rather than one reconstructed by a few operations. Which floor applies
// is a property of the quantity, not of the assertion's form: REAL_REL_FLOOR is a few epsilon and
// is unreachable for anything the solver integrated. In a double build the floor is zero, so this
// reduces exactly to EXPECT_REAL_REL.
#define EXPECT_REAL_SOLVE_REL(actual, reference, rel_tol) \
  EXPECT_NEAR((actual), (reference), micm::test::SolveRelativeTolerance(rel_tol) * std::abs((double)(reference)))

// Agreement between two evaluations of the same closed-form expression.
//   double build -> EXPECT_DOUBLE_EQ, so double-precision tests keep the tolerance they had
//   float  build -> a small multiple of epsilon, relative to the reference
//
// Use in place of EXPECT_REAL_EQ when the two sides compute the same formula by different routes --
// a test re-deriving a rate constant against micm's own implementation, for instance. Not for
// values legitimately near zero.
#ifdef MICM_USE_SINGLE
  #define EXPECT_REAL_FORMULA_EQ(a, b) EXPECT_NEAR((a), (b), micm::test::FormulaTolerance((double)(b)))
#else
  #define EXPECT_REAL_FORMULA_EQ(a, b) EXPECT_DOUBLE_EQ((a), (b))
#endif
