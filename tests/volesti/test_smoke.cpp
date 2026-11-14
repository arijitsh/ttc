// Smoke tests for the scaffold + the user-reported long-box failure.
#include "doctest.h"
#include "common.hpp"
#include <iostream>

using namespace vt;

static double runTtc(const Body& B)
{
  HPoly P = toHPoly(B);
  volestiBall(P);
  return ttcVolume(P);
}

TEST_CASE("smoke: unit cube dims 2..6 within eps/12-ish")
{
  for (int n = 2; n <= 6; ++n)
  {
    Body B = cube(n);
    double v = runTtc(B);
    std::cout << B.name << " est=" << v << " exact=" << B.exactVolume
              << " rel=" << v / B.exactVolume - 1 << "\n";
    CHECK(std::abs(v / B.exactVolume - 1) < 0.2);
  }
}

TEST_CASE("user report: long box [0,0]-[2,1000]")
{
  Body B = box({0, 0}, {2, 1000});
  double v = runTtc(B);
  std::cout << "long box est=" << v << " exact=" << B.exactVolume
            << " rel=" << v / B.exactVolume - 1 << "\n";
  CHECK(std::abs(v / B.exactVolume - 1) < 0.2);
}
