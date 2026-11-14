// Accuracy statistics of the volesti Gaussian-cooling volume oracle on a
// catalogue of standard, mostly well-rounded bodies with exact volumes.
//
// For each body we run
//   * ttcVolume(P)                         -- ttc's exact deterministic call (rounding + oracle)
//                                             (epsilon/12 = 0.0667, walk 10+n^2,
//                                              VolumeRNG seeded with 3)
//   * gaussianOracle(P, seed, 0.8/12, 10+n^2) over >= 20 seeds
// and report mean relative bias, sd, and the fraction of runs within the
// oracle target (+-0.0667) and twice the target (+-0.1333).
//
// Assertions (they must hold for a correct oracle; failures expose bugs):
//   * |ttcVolume/exact - 1|                  <= 2 * 0.0667
//   * |mean_seeds(est/exact - 1)|            <  0.0667
//   * fraction of seeds with |rel| <= 0.1333 >= 0.80
//
// The inner ball is pre-set (exact Chebyshev ball where known, otherwise
// volesti's own LP ball), mirroring ttc's reliableInnerBall.
//
// Fast cases run by default (< ~2 min).  Large dimensions are behind the
// test-suite "slow" (skipped by default; run with --no-skip or -ts=slow).
#include "doctest.h"
#include "common.hpp"

#include <chrono>
#include <cstdio>
#include <iostream>

using namespace vt;

namespace
{

constexpr double kTarget = 0.8 / 12.0;  // ttc's oracle error (epsilon = 0.8)
constexpr double kTol2 = 2.0 * kTarget;

// ---------------------------------------------------------------- bodies
struct Case
{
  Body body;
  bool exactBall = false;  // use (center, radius) below, else volesti's LP ball
  VT center;
  double radius = 0.0;
  bool normalizeRows = true;  // false: keep the user's row scaling
};

Case withBall(Body B, const VT& c, double r)
{
  Case C;
  C.body = std::move(B);
  C.exactBall = true;
  C.center = c;
  C.radius = r;
  return C;
}

Case lpBall(Body B)
{
  Case C;
  C.body = std::move(B);
  return C;
}

HPoly prepare(const Case& C)
{
  HPoly P = toHPoly(C.body);
  if (C.exactBall)
  {
    if (C.normalizeRows) P.normalize();
    P.set_InnerBall({Point(C.center), C.radius});
  }
  else
  {
    volestiBall(P);
  }
  return P;
}

Case cubeCase(int n) { return withBall(cube(n), VT::Zero(n), 1.0); }

// Standard simplex: Chebyshev ball centre r*1, r = 1/(n + sqrt n).
Case simplexCase(int n)
{
  const double r = 1.0 / (n + std::sqrt(double(n)));
  return withBall(simplex(n), VT::Constant(n, r), r);
}

Case crossCase(int n) { return withBall(crossPolytope(n), VT::Zero(n), 1.0 / std::sqrt(double(n))); }

Case rotatedCubeCase(int n, unsigned seed)
{
  Body B = linearImage(cube(n), randomRotation(n, seed), VT::Zero(n), "rot");
  return withBall(B, VT::Zero(n), 1.0);
}

// Delta_k x Delta_l (standard simplices), volume 1/(k! l!).
Case simplexProduct(int k, int l)
{
  const int n = k + l;
  Body S1 = simplex(k), S2 = simplex(l);
  Body B;
  B.name = "simplex" + std::to_string(k) + "xsimplex" + std::to_string(l);
  B.A = MT::Zero(S1.A.rows() + S2.A.rows(), n);
  B.b = VT::Zero(B.A.rows());
  B.A.block(0, 0, S1.A.rows(), k) = S1.A;
  B.A.block(S1.A.rows(), k, S2.A.rows(), l) = S2.A;
  B.b << S1.b, S2.b;
  B.exactVolume = S1.exactVolume * S2.exactVolume;
  const double r1 = 1.0 / (k + std::sqrt(double(k))), r2 = 1.0 / (l + std::sqrt(double(l)));
  VT c(n);
  c << VT::Constant(k, r1), VT::Constant(l, r2);
  return withBall(B, c, std::min(r1, r2));
}

// Regular m-gon with inradius 1: area m tan(pi/m).
Case regularPolygon(int m)
{
  Body B;
  B.name = "regular" + std::to_string(m) + "gon";
  B.A = MT(m, 2);
  B.b = VT::Ones(m);
  for (int j = 0; j < m; ++j)
  {
    const double t = 2.0 * M_PI * j / m;
    B.A(j, 0) = std::cos(t);
    B.A(j, 1) = std::sin(t);
  }
  B.exactVolume = m * std::tan(M_PI / m);
  return withBall(B, VT::Zero(2), 1.0);
}

// Birkhoff polytope B_3 in the coordinates (x11, x12, x21, x22) of its free
// 2x2 block.  Its volume there is 1/8 (the "9/8" relative volume of B_3 is
// w.r.t. the lattice of the affine hull, index 3^2 = 9 bigger).  1/8 was
// cross-checked with 2e8 uniform samples in [0,1]^4: 0.125040.
Case birkhoff3()
{
  Body B;
  B.name = "birkhoff3";
  B.A = MT::Zero(9, 4);
  B.b = VT::Zero(9);
  // x_ij >= 0
  for (int i = 0; i < 4; ++i) B.A(i, i) = -1.0;
  // x13 = 1 - x11 - x12 >= 0, x23 = 1 - x21 - x22 >= 0
  B.A.row(4) << 1, 1, 0, 0;  B.b(4) = 1;
  B.A.row(5) << 0, 0, 1, 1;  B.b(5) = 1;
  // x31 = 1 - x11 - x21 >= 0, x32 = 1 - x12 - x22 >= 0
  B.A.row(6) << 1, 0, 1, 0;  B.b(6) = 1;
  B.A.row(7) << 0, 1, 0, 1;  B.b(7) = 1;
  // x33 = x11 + x12 + x21 + x22 - 1 >= 0
  B.A.row(8) << -1, -1, -1, -1;  B.b(8) = -1;
  B.exactVolume = 1.0 / 8.0;
  return lpBall(B);
}

// Box [0,1]x[0,2]x...x[0,n] with redundant rows appended: duplicated facets,
// facets far away, a facet scaled-duplicate and a redundant facet touching
// a vertex.
Case redundantBox(int n)
{
  std::vector<double> lo(n, 0.0), hi(n);
  for (int i = 0; i < n; ++i) hi[i] = 1.0 + i;
  Body B0 = box(lo, hi);
  const int extra = 2 * n + 4;
  Body B;
  B.name = "redundantbox" + std::to_string(n);
  B.A = MT::Zero(B0.A.rows() + extra, n);
  B.b = VT::Zero(B.A.rows());
  B.A.topRows(B0.A.rows()) = B0.A;
  B.b.head(B0.A.rows()) = B0.b;
  int r = B0.A.rows();
  for (int i = 0; i < n; ++i)  // exact duplicates of the upper facets
  {
    B.A.row(r) = B0.A.row(2 * i);
    B.b(r++) = B0.b(2 * i);
  }
  for (int i = 0; i < n; ++i)  // scaled duplicates of the lower facets (x 7)
  {
    B.A.row(r) = 7.0 * B0.A.row(2 * i + 1);
    B.b(r++) = 7.0 * B0.b(2 * i + 1);
  }
  B.A.row(r).setOnes();  // far: sum x <= 1e3
  B.b(r++) = 1e3;
  B.A.row(r).setZero();  // far: x_0 <= 1e4
  B.A(r, 0) = 1.0;
  B.b(r++) = 1e4;
  B.A.row(r).setOnes();  // touching the top vertex: sum x <= sum hi
  B.b(r++) = std::accumulate(hi.begin(), hi.end(), 0.0);
  B.A.row(r).setZero();  // far, negative side: -x_{n-1} <= 50
  B.A(r, n - 1) = -1.0;
  B.b(r++) = 50.0;
  B.exactVolume = B0.exactVolume;
  return lpBall(B);
}

// Cube [-1,1]^n given with rows multiplied alternately by 1e6 and 1e-6.
// The ball is set WITHOUT normalizing the rows first, i.e. the polytope
// reaches the oracle with the user's (unnormalized) scaling.
Case rowScaledCube(int n)
{
  Body B = cube(n);
  B.name = "rowscaledcube" + std::to_string(n);
  for (int i = 0; i < B.A.rows(); ++i)
  {
    const double s = (i % 2 == 0) ? 1e6 : 1e-6;
    B.A.row(i) *= s;
    B.b(i) *= s;
  }
  Case C = withBall(B, VT::Zero(n), 1.0);
  C.normalizeRows = false;
  return C;
}

// Simplex given with rows scaled by 1e6 (facets x_i>=0) and 1e-6 (sum<=1),
// inner ball from volesti's LP.
Case rowScaledSimplex(int n)
{
  Body B = simplex(n);
  B.name = "rowscaledsimplex" + std::to_string(n);
  for (int i = 0; i < n; ++i) B.A.row(i) *= 1e6;
  B.A.row(n) *= 1e-6;
  B.b(n) *= 1e-6;
  return lpBall(B);
}

// ------------------------------------------------------------ measurement
struct Result
{
  double ttc = 0, ttcRel = 0;
  Stats S, S2;  // within target, within 2x target
};

Result measure(const Case& C, int seeds = 20)
{
  const HPoly P = prepare(C);
  const int n = static_cast<int>(P.dimension());
  const unsigned wl = ttcWalkLength(n);
  Result R;
  R.ttc = ttcVolume(P);  // ttc's full pipeline: rounding + oracle
  R.ttcRel = R.ttc / C.body.exactVolume - 1.0;

  // Seeds run sequentially (one compute thread: the machine is shared).
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<double> est(seeds, 0.0);
  for (int s = 0; s < seeds; ++s) est[s] = gaussianOracle(P, 1000u + 7919u * s, kTarget, wl);
  const double secs =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  R.S = relStats(est, C.body.exactVolume, kTarget);
  R.S2 = relStats(est, C.body.exactVolume, kTol2);
  char buf[512];
  std::snprintf(buf, sizeof buf,
                "%-26s n=%2d exact=%-11.5g ttc=%-11.5g ttcRel=%+.4f | %s within2x=%.2f (%.1fs)",
                C.body.name.c_str(), n, C.body.exactVolume, R.ttc, R.ttcRel, fmt(R.S).c_str(),
                R.S2.fracWithin, secs);
  std::cout << buf << std::endl;
  return R;
}

void checkCase(const Case& C, int seeds = 20)
{
  Result R = measure(C, seeds);
  INFO(C.body.name);
  CHECK(std::isfinite(R.ttc));
  CHECK(std::abs(R.ttcRel) <= kTol2);
  CHECK(R.S.nonFinite == 0);
  CHECK(std::abs(R.S.meanRel) < kTarget);
  CHECK(R.S2.fracWithin >= 0.8);
}

}  // namespace

// ================================================================== fast
// Budget: one thread, < ~2 min in total.  Everything else is in "slow".
TEST_CASE("shapes: interval (dimension 1)")
{
  // ttc can produce 1-dimensional polytopes.  With d = 1 the cooling ratio
  // 1 - 1/d is 0, so the schedule jumps straight from a_0 to the uniform
  // (a_vals = {a_0, 0}).  Measured: works, bias < 1%, sd ~ 3%.
  checkCase(cubeCase(1));
  checkCase(withBall(box({0.0}, {5.0}), VT::Constant(1, 2.5), 2.5));
}

TEST_CASE("shapes: small cubes, cross-polytopes, rotated cubes")
{
  for (int n = 2; n <= 5; ++n) checkCase(cubeCase(n));
  checkCase(crossCase(3));
  checkCase(rotatedCubeCase(3, 17u + 3));
}

TEST_CASE("shapes: standard simplices n=2..6")
{
  for (int n = 2; n <= 6; ++n) checkCase(simplexCase(n));
}

TEST_CASE("shapes: regular polygons")
{
  checkCase(regularPolygon(6));
  checkCase(regularPolygon(8));
  checkCase(regularPolygon(3));
}

TEST_CASE("shapes: Birkhoff B_3")
{
  checkCase(birkhoff3());
}

TEST_CASE("shapes: boxes with redundant facets n=2,3")
{
  for (int n : {2, 3}) checkCase(redundantBox(n));
}

TEST_CASE("shapes: unnormalized row scales (1e6 / 1e-6)")
{
  checkCase(rowScaledCube(2));
  checkCase(rowScaledSimplex(3));
}

// ================================================================== slow
// Run with:  bin/test_standard_shapes -ts=slow --no-skip
TEST_SUITE("slow")
{
  TEST_CASE("shapes: cubes n=6..20" * doctest::skip())
  {
    for (int n = 6; n <= 20; ++n) checkCase(cubeCase(n));
  }

  // R/r of the standard simplex grows like n; the spread of the estimate
  // grows with it (sd 0.047 at n=2 -> 0.10 at n=7) and from n~6 the
  // +-2x-target coverage drops below 80% (FAILS at n=7,8 on HEAD).
  TEST_CASE("shapes: standard simplices n=7..15" * doctest::skip())
  {
    for (int n = 7; n <= 15; ++n) checkCase(simplexCase(n));
  }

  TEST_CASE("shapes: cross-polytopes n=2..8" * doctest::skip())
  {
    for (int n = 2; n <= 8; ++n) checkCase(crossCase(n));
  }

  TEST_CASE("shapes: rotated cubes n=2..10" * doctest::skip())
  {
    for (int n : {2, 4, 5, 6, 8, 10}) checkCase(rotatedCubeCase(n, 17u + n));
  }

  TEST_CASE("shapes: products of simplices" * doctest::skip())
  {
    checkCase(simplexProduct(2, 2));
    checkCase(simplexProduct(3, 3));
    checkCase(simplexProduct(2, 4));
    checkCase(simplexProduct(4, 4));
    checkCase(simplexProduct(5, 5));
  }

  TEST_CASE("shapes: boxes with redundant facets n=5" * doctest::skip())
  {
    checkCase(redundantBox(5));
  }

  TEST_CASE("shapes: unnormalized row scales n=4" * doctest::skip())
  {
    checkCase(rowScaledCube(4));
  }
}
