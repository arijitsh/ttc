// Elongated / anisotropic bodies: failure envelope of ttc's volume oracle.
//
// ttc calls volume_cooling_gaussians<GaussianBallWalk>(P, eps/12, 10+n^2) on
// the *unrounded* polytope.  The Gaussian ball walk's step is
//     delta = 4 r / sqrt(max(1,a) n)          (gaussian_ball_walk.hpp:52,
//                                              volume_cooling_gaussians.hpp:238,396)
// where r is the inner radius.  In a body of diameter D >> r the walk needs
// ~(D/delta)^2 steps to spread, but gets 10+n^2, so every annealing ratio is
// estimated from points stuck near the centre -> large underestimates.  The
// formula also mixes units (r is a length, a is 1/length^2), so a uniformly
// *small* body is mis-estimated too (thin axis / scaled cube cases).
//
// All "bug" tests below assert correct behaviour (|rel err| within a few times
// the 0.0667 oracle target) and FAIL on the current oracle; they pass once the
// oracle rounds the body first (see the report / validated overlay patch).
// "control" tests pass today and must keep passing.
// Slow sweeps are behind doctest::skip(); run with  -ts=slow  --no-skip.
#include "doctest.h"
#include "common.hpp"
#include <cstdio>
#include <iostream>

using namespace vt;

namespace
{

// Single-run tolerance: 3x the oracle target eps/12 = 0.0667 (unit cubes in
// dims 2..10 have per-run sd 0.03-0.06 on the current oracle).
constexpr double kTol = 0.2;
constexpr double kTarget = 0.8 / 12.0;

double runTtc(const Body& B)
{
  HPoly P = toHPoly(B);
  volestiBall(P);
  return ttcVolume(P);
}

double relErr(double est, const Body& B) { return est / B.exactVolume - 1.0; }

// [0,1]^n with axis n-1 stretched to [0,R].
Body longBox(int n, double R)
{
  std::vector<double> lo(n, 0.0), hi(n, 1.0);
  hi[n - 1] = R;
  Body B = box(lo, hi);
  B.name = "longbox" + std::to_string(n) + "_R" + std::to_string((long long)R);
  return B;
}

// [0,1]^n with the last n/2 axes stretched to [0,R].
Body manyLongBox(int n, double R)
{
  std::vector<double> lo(n, 0.0), hi(n, 1.0);
  for (int i = n / 2; i < n; ++i) hi[i] = R;
  Body B = box(lo, hi);
  B.name = "manylong" + std::to_string(n) + "_R" + std::to_string((long long)R);
  return B;
}

// [0,w] x [0,1]^{n-1}.
Body thinBox(int n, double w)
{
  std::vector<double> lo(n, 0.0), hi(n, 1.0);
  hi[0] = w;
  Body B = box(lo, hi);
  B.name = "thin" + std::to_string(n) + "_w" + std::to_string(w);
  return B;
}

// Standard simplex with the last coordinate stretched by R.
Body longSimplex(int n, double R)
{
  MT Q = MT::Identity(n, n);
  Q(n - 1, n - 1) = R;
  return linearImage(simplex(n), Q, VT::Zero(n), "stretch" + std::to_string((long long)R));
}

Body rotated(const Body& B, unsigned seed)
{
  const int n = static_cast<int>(B.A.cols());
  return linearImage(B, randomRotation(n, seed), VT::Zero(n), "rot" + std::to_string(seed));
}

void report(const Body& B, double est)
{
  std::printf("  %-34s est=%-12.6g exact=%-12.6g rel=%+.4f\n", B.name.c_str(), est,
              B.exactVolume, relErr(est, B));
}

Stats seeds(const Body& B, int k, unsigned wl, double tol = kTol)
{
  std::vector<double> est;
  for (int s = 0; s < k; ++s)
  {
    HPoly P = toHPoly(B);
    volestiBall(P);
    est.push_back(gaussianOracle(P, 1000u + s, kTarget, wl));
  }
  Stats S = relStats(est, B.exactVolume, tol);
  std::printf("  %-34s wl=%-5u %s\n", B.name.c_str(), wl, fmt(S).c_str());
  return S;
}

// Same statistics through ttc's full pipeline (rounding + oracle).
Stats seedsPipeline(const Body& B, int k, double tol = kTol)
{
  std::vector<double> est;
  for (int s = 0; s < k; ++s)
  {
    HPoly P = toHPoly(B);
    volestiBall(P);
    est.push_back(ttcVolumeSeeded(P, 1000u + s));
  }
  Stats S = relStats(est, B.exactVolume, tol);
  std::printf("  %-34s pipeline %s\n", B.name.c_str(), fmt(S).c_str());
  return S;
}

}  // namespace

// ------------------------------------------------------------------ controls

TEST_CASE("elongated control: aspect ratio 1-2 boxes are accurate (n=2,3,5)")
{
  for (int n : {2, 3, 5})
    for (double R : {1.0, 2.0})
    {
      Body B = longBox(n, R);
      double v = runTtc(B);
      report(B, v);
      CHECK(std::abs(relErr(v, B)) < kTol);
    }
}

TEST_CASE("elongated control: the long box mapped back to a cube is accurate")
{
  // Same body as the [0,1]x[0,1000] failure, but the caller applies the
  // diagonal map first and multiplies by det: this is what rounding buys.
  Body B = longBox(2, 1000.0);
  MT D = MT::Identity(2, 2);
  D(1, 1) = 1.0 / 1000.0;
  Body C = linearImage(B, D, VT::Zero(2), "unstretched");
  double v = runTtc(C) * 1000.0;
  report(B, v);
  CHECK(std::abs(v / B.exactVolume - 1.0) < kTol);
}

// --------------------------------------------------------- bug: long boxes

TEST_CASE("elongated: 2-D long box [0,1]x[0,R], R = 10..1e4")
{
  // Current oracle: R=10 -0.10, R=100 -0.67, R=1000 -0.94, R=1e4 -0.997
  // (mean over 10 seeds).  Error grows monotonically, always downwards.
  for (double R : {10.0, 100.0, 1000.0, 10000.0})
  {
    Body B = longBox(2, R);
    double v = runTtc(B);
    report(B, v);
    CHECK_MESSAGE(std::abs(relErr(v, B)) < kTol, B.name);
  }
}

TEST_CASE("elongated: one long axis, n=3 and n=5, R=100")
{
  // Current oracle (10 seeds): n=3 -0.73, n=5 -0.81.
  for (int n : {3, 5})
  {
    Body B = longBox(n, 100.0);
    double v = runTtc(B);
    report(B, v);
    CHECK_MESSAGE(std::abs(relErr(v, B)) < kTol, B.name);
  }
}

TEST_CASE("elongated: several long axes (n=5, axes 2..4 of length 300)")
{
  // Current oracle (10 seeds): -0.993; every axis long compounds the error.
  Body B = manyLongBox(5, 300.0);
  double v = runTtc(B);
  report(B, v);
  CHECK(std::abs(relErr(v, B)) < kTol);
}

TEST_CASE("elongated: statistics across seeds, 3-D long box R=30")
{
  // Moderate aspect: the mean is biased low AND the spread is huge
  // (current: meanRel -0.57, sd 0.18, 0/10 within 0.0667).
  Body B = longBox(3, 30.0);
  Stats S = seedsPipeline(B, 6);
  CHECK(std::abs(S.meanRel) < 0.1);
  CHECK(S.fracWithin >= 0.8);
}

// ------------------------------------------------- bug: thin axis / scaling

TEST_CASE("elongated: one thin axis [0,1e-3]x[0,1]^{n-1}, n=2,3,5")
{
  // Current oracle: rel = -0.9998 (n=2), -1.0000 (n=3,5): returns ~0.
  // Same shape up to scale as a long box, but worse: delta=4r/sqrt(a n)
  // is ~sigma*r, so for r=5e-4 the walk is frozen and the annealing
  // jumps a_0 ~ 1.5e5 -> a_1 ~ 0.35 in a single step.
  for (int n : {2, 3, 5})
  {
    Body B = thinBox(n, 1e-3);
    double v = runTtc(B);
    report(B, v);
    CHECK_MESSAGE(std::abs(relErr(v, B)) < kTol, B.name);
  }
}

TEST_CASE("elongated: scale invariance, cube [0,s]^3 for s = 1e-3, 1e-2")
{
  // An isotropic body that is merely small: current rel -0.93 (s=1e-3),
  // -0.82 (s=1e-2) vs +0.01 for s=1.  The oracle is not scale invariant.
  for (double s : {1e-3, 1e-2})
  {
    Body B = cube(3, 0.0, s);
    double v = runTtc(B);
    report(B, v);
    CHECK_MESSAGE(std::abs(relErr(v, B)) < kTol, B.name);
  }
}

// --------------------------------------- bug: simplices and rotated bodies

TEST_CASE("elongated: stretched simplex n=3, factor 100")
{
  // Current oracle (10 seeds): meanRel -0.60, sd 0.33.
  Body B = longSimplex(3, 100.0);
  double v = runTtc(B);
  report(B, v);
  CHECK(std::abs(relErr(v, B)) < kTol);
}

TEST_CASE("elongated: rotated long box n=2 R=1000 and n=5 R=100")
{
  // Rotation removes any axis alignment (rules out coordinate-walk fixes).
  // Current (10 seeds): n=2 R=1000 -0.96.
  for (auto B : {rotated(longBox(2, 1000.0), 7), rotated(longBox(5, 100.0), 7)})
  {
    double v = runTtc(B);
    report(B, v);
    CHECK_MESSAGE(std::abs(relErr(v, B)) < kTol, B.name);
  }
}

// --------------------------------------------------------------- slow sweeps

TEST_CASE("elongated sweep: walk length 1x/10x/100x on 2-D box R=100" * doctest::test_suite("slow")
          * doctest::skip())
{
  // Distinguishes mixing from structure.  Measured (4 seeds):
  //   wl=14   meanRel -0.63 sd 0.31   (2.7 s)
  //   wl=140  meanRel -0.20 sd 0.25   (39 s)
  //   wl=1400 meanRel -0.08 sd 0.18   (513 s)
  // => error shrinks with walk length (a mixing problem), but the needed
  // length grows like (R/r)^2, so no fixed 10+n^2 walk length can work.
  Body B = longBox(2, 100.0);
  double prev = -1e9;
  for (unsigned mult : {1u, 10u, 100u})
  {
    Stats S = seeds(B, 4, mult * ttcWalkLength(2));
    CHECK(S.meanRel > prev - 0.05);  // non-decreasing (bias shrinks)
    prev = S.meanRel;
  }
  CHECK(std::abs(prev) < 0.15);
}

TEST_CASE("elongated sweep: aspect ratio x dimension table" * doctest::test_suite("slow")
          * doctest::skip())
{
  // Prints meanRel over 5 seeds; asserts every cell is within 2x target.
  for (int n : {2, 3, 5, 8, 10})
    for (double R : {1.0, 2.0, 5.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 1e4})
    {
      for (const Body& B : {longBox(n, R), manyLongBox(n, R), thinBox(n, 1.0 / R),
                            longSimplex(n, R)})
      {
        Stats S = seeds(B, 5, ttcWalkLength(n));
        CHECK_MESSAGE(std::abs(S.meanRel) < 2 * kTarget, B.name);
      }
    }
}
