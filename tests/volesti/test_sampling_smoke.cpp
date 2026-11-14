// Smoke test for ttc's single-polytope sampler: are the samples spread over
// the whole body?  For a box the uniform distribution makes each coordinate
// uniform on its side, so the fraction of samples in any slab is known.
#include "doctest.h"
#include "common.hpp"
#include <cstdio>

using namespace vt;

// Fraction of points whose coordinate j lies in [lo, hi].
static double frac(const std::vector<VT>& pts, int j, double lo, double hi)
{
  int k = 0;
  for (const auto& p : pts) k += (p(j) >= lo && p(j) <= hi);
  return double(k) / pts.size();
}

TEST_CASE("sampling smoke: unit square, billiard proposals cover all four quadrants")
{
  Body B = box({0, 0}, {1, 1});
  auto pts = ttcProposals(B.A, B.b, 4000, 7);
  for (double lo : {0.0, 0.5})
  {
    const double f = frac(pts, 0, lo, lo + 0.5);
    std::printf("unit square x in [%.1f,%.1f]: %.3f (uniform 0.5)\n", lo, lo + 0.5, f);
    CHECK(std::abs(f - 0.5) < 0.05);
  }
}

TEST_CASE("sampling smoke: long box [0,2]x[0,1000], samples must reach both ends")
{
  // Uniform: 10% of the samples in each tenth of the long axis.
  Body B = box({0, 0}, {2, 1000});
  auto pts = ttcProposals(B.A, B.b, 4000, 7);
  for (double lo : {0.0, 450.0, 900.0})
  {
    const double f = frac(pts, 1, lo, lo + 100);
    std::printf("long box y in [%g,%g]: %.3f (uniform 0.100)\n", lo, lo + 100, f);
    CHECK(std::abs(f - 0.1) < 0.03);
  }
}
