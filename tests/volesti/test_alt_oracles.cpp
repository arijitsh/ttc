// Compare volesti's alternative volume oracles on the same bodies.
//
//   CG-Ball   volume_cooling_gaussians<GaussianBallWalk>  err=eps/12, wl=10+n^2  (what ttc uses)
//   CG-CDHR   volume_cooling_gaussians<GaussianCDHRWalk>  err=eps/12, wl=10+n/10 (volesti default walk)
//   CB-Bill   volume_cooling_balls<BilliardWalk>          err=eps/12, wl=1
//   CB-ABill  volume_cooling_balls<AcceleratedBilliardWalk> err=eps/12, wl=1
//   CB-CDHR   volume_cooling_balls<CDHRWalk>              err=eps/12, wl=10+n/10
//   SOB-CDHR / SOB-RDHR / SOB-Bill   volume_sequence_of_balls  err=0.2, wl=10+n/10
//   R+<x>     max-inscribed-ellipsoid rounding (inscribed_ellipsoid_rounding; d=2 falls back
//             to min_sampling_covering_ellipsoid_rounding) first,
//             then oracle <x> on the rounded body, times |det T|.
//
// eps = 0.8 (ttc default) => oracle error target 0.0667.
//
// Test cases:
//   * "alt oracles: ..." (fast) assert accuracy of the oracles that should be accurate, over
//     several seeds, on the skinny bodies that break ttc's oracle.
//   * "alt oracles: full sweep" (tagged slow, skipped by default) prints the full table
//     (>= 10 seeds per body/oracle, runtime).  Run with:
//        bin/test_alt_oracles -ts=slow -nskip
//     or a subset (comma-separated name prefixes):
//        ALT_ORACLES=CG-Ball,R+CB ALT_BODIES=box2-,cube5 ALT_SEEDS=10 ALT_BUDGET=60 ...
#include "doctest.h"
#include "common.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>

#include "volume/volume_cooling_balls.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"
#include "preprocess/min_sampling_covering_ellipsoid_rounding.hpp"

using namespace vt;

namespace
{

constexpr double kEps = 0.8;
constexpr double kErr = kEps / 12.0;  // 0.0667, the oracle target ttc requests

enum class Oracle
{
  CG_Ball,
  CG_CDHR,
  CB_Bill,
  CB_ABill,
  CB_CDHR,
  SOB_CDHR,
  SOB_RDHR,
  SOB_Bill,
};

const std::vector<std::pair<Oracle, std::string>> kOracles = {
    {Oracle::CG_Ball, "CG-Ball"},   {Oracle::CG_CDHR, "CG-CDHR"},   {Oracle::CB_Bill, "CB-Bill"},
    {Oracle::CB_ABill, "CB-ABill"}, {Oracle::CB_CDHR, "CB-CDHR"},   {Oracle::SOB_CDHR, "SOB-CDHR"},
    {Oracle::SOB_RDHR, "SOB-RDHR"}, {Oracle::SOB_Bill, "SOB-Bill"},
};

std::string oracleName(Oracle o)
{
  for (auto& p : kOracles)
    if (p.first == o) return p.second;
  return "?";
}

using RNG = TtcSampleRNG;  // BoostRandomNumberGenerator<boost::mt19937,double>, seedable

double callOracle(Oracle o, HPoly& P, unsigned seed)
{
  const unsigned n = P.dimension();
  RNG rng(static_cast<int>(n));
  rng.set_seed(seed);
  const unsigned wlCdhr = 10 + n / 10;
  switch (o)
  {
    case Oracle::CG_Ball:
      return volume_cooling_gaussians<GaussianBallWalk>(P, rng, kErr, ttcWalkLength(n));
    case Oracle::CG_CDHR:
      return volume_cooling_gaussians<GaussianCDHRWalk>(P, rng, kErr, wlCdhr);
    case Oracle::CB_Bill:
      return volume_cooling_balls<BilliardWalk>(P, rng, kErr, 1u).second;
    case Oracle::CB_ABill:
      return volume_cooling_balls<AcceleratedBilliardWalk>(P, rng, kErr, 1u).second;
    case Oracle::CB_CDHR:
      return volume_cooling_balls<CDHRWalk>(P, rng, kErr, wlCdhr).second;
    // SOB with err=0.0667 needs 400 n ln n / err^2 points per ball (~10^6 at n=5): use 0.2.
    case Oracle::SOB_CDHR:
      return volume_sequence_of_balls<CDHRWalk>(P, rng, 0.2, wlCdhr);
    case Oracle::SOB_RDHR:
      return volume_sequence_of_balls<RDHRWalk>(P, rng, 0.2, wlCdhr);
    case Oracle::SOB_Bill:
      return volume_sequence_of_balls<BilliardWalk>(P, rng, 0.2, 1u);
  }
  return -1;
}

struct Run
{
  double vol;
  double secs;
};

// Rounding: volesti's max-inscribed-ellipsoid rounding (inscribed_ellipsoid_rounding,
// deterministic).  It throws for d=2: max_inscribed_ellipsoid builds a Spectra
// SymEigsSolver with nev=2, ncv=min(max(10,d/5),d) (rounding_util_functions.hpp:245,
// max_inscribed_ellipsoid.hpp:132), and Spectra requires nev <= d-1.  In that case fall
// back to min_sampling_covering_ellipsoid_rounding<CDHRWalk> (Khachiyan on CDHR samples).
// Returns |det T|: volume(original) = |det T| * volume(rounded P).
double roundPolytope(HPoly& P, unsigned seed)
{
  auto ball = P.ComputeInnerBall();
  HPoly backup(P);
  try
  {
    return std::get<2>(inscribed_ellipsoid_rounding<MT, VT, NT>(P, ball.first));
  }
  catch (const std::exception&)
  {
    P = backup;
    RNG rng(static_cast<int>(P.dimension()));
    rng.set_seed(seed);
    return std::get<2>(min_sampling_covering_ellipsoid_rounding<CDHRWalk, MT, VT>(
        P, ball, 10 + 10 * P.dimension(), rng));
  }
}

// Run oracle o on body B (optionally rounded first).
Run runOracle(Oracle o, const Body& B, unsigned seed, bool round)
{
  auto t0 = std::chrono::steady_clock::now();
  HPoly P = toHPoly(B);
  double scale = 1.0;
  double v;
  if (round)
  {
    scale = roundPolytope(P, seed);
    HPoly Q(P.dimension(), P.get_mat(), P.get_vec());
    Q.ComputeInnerBall();
    v = callOracle(o, Q, seed);
  }
  else
  {
    P.ComputeInnerBall();  // ttc pre-sets the (exact) Chebyshev ball; same here
    v = callOracle(o, P, seed);
  }
  if (v > 0 && std::isfinite(v)) v *= scale;
  double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return {v, secs};
}

// One long side: [0,1]^(n-1) x [0,a]  (volume a).
Body longBox(int n, double aspect)
{
  std::vector<double> lo(n, 0.0), hi(n, 1.0);
  hi[n - 1] = aspect;
  Body B = box(lo, hi);
  B.name = "box" + std::to_string(n) + "-asp" + std::to_string(int(aspect));
  return B;
}

Body userBox()
{
  Body B = box({0, 0}, {2, 1000});
  B.name = "box[0,2]x[0,1000]";
  return B;
}

Body rotatedLongBox5()
{
  Body B = longBox(5, 100.0);
  VT t = VT::Constant(5, 0.5);
  Body R = linearImage(B, randomRotation(5, 7), t, "rot");
  R.name = "rot-box5-asp100";
  return R;
}

std::vector<Body> allBodies()
{
  std::vector<Body> v;
  v.push_back(userBox());
  for (int n : {2, 5, 10})
    for (double a : {10.0, 100.0, 1000.0}) v.push_back(longBox(n, a));
  for (int n : {5, 10, 15, 20})
  {
    Body c = cube(n);
    c.name = "cube" + std::to_string(n);
    v.push_back(c);
  }
  v.push_back(simplex(5));
  v.push_back(simplex(10));
  v.push_back(rotatedLongBox5());
  return v;
}

struct Cell
{
  Stats S;
  double meanSecs = 0;
  int runs = 0;
};

Cell evalCell(Oracle o, const Body& B, bool round, int seeds, double tol,
              double timeBudget = 1e9)
{
  std::vector<double> est;
  double tot = 0;
  for (int s = 1; s <= seeds; ++s)
  {
    Run r = runOracle(o, B, 1000u + s, round);
    est.push_back(r.vol);
    tot += r.secs;
    if (tot > timeBudget) break;  // keep the slow sweep bounded
  }
  Cell c;
  c.S = relStats(est, B.exactVolume, tol);
  c.runs = static_cast<int>(est.size());
  c.meanSecs = tot / c.runs;
  return c;
}

std::string label(Oracle o, bool round) { return (round ? "R+" : "") + oracleName(o); }

void report(const Body& B, Oracle o, bool round, const Cell& c)
{
  std::printf("%-20s %-12s %s  t=%.3fs\n", B.name.c_str(), label(o, round).c_str(),
              fmt(c.S).c_str(), c.meanSecs);
  std::fflush(stdout);
}

// Accuracy bound used in assertions: oracle target is 0.0667; we allow 3x on the
// mean relative error across seeds and require >= 80% of runs within 3x target.
constexpr double kTol = 3 * kErr;  // 0.2

void checkAccurate(const Body& B, Oracle o, bool round, int seeds)
{
  Cell c = evalCell(o, B, round, seeds, kTol);
  report(B, o, round, c);
  INFO(B.name << " " << label(o, round) << " " << fmt(c.S));
  CHECK(c.S.nonFinite == 0);
  CHECK(std::abs(c.S.meanRel) <= kTol);
  CHECK(c.S.fracWithin >= 0.8);
}

bool selected(const char* env, const std::string& name)
{
  const char* s = std::getenv(env);
  if (!s || !*s) return true;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ','))
    if (!tok.empty() && name.rfind(tok, 0) == 0) return true;
  return false;
}

}  // namespace

// ---------------------------------------------------------------------------- fast

TEST_CASE("alt oracles: ttc's CG-Ball on skinny / rotated bodies is accurate")
{
  // Bug reproducer: user box exact 2000, ttc config returns ~7.8 (meanRel -0.87 over 10
  // seeds); rot-box5-asp100 meanRel -0.72.  Must pass once the oracle is fixed.
  checkAccurate(userBox(), Oracle::CG_Ball, false, 5);
  checkAccurate(longBox(2, 100), Oracle::CG_Ball, false, 3);
  checkAccurate(rotatedLongBox5(), Oracle::CG_Ball, false, 3);
}

TEST_CASE("alt oracles: rounding + {CG-Ball, CG-CDHR, CB-Bill} on skinny/rotated bodies")
{
  // Max-inscribed-ellipsoid rounding makes every oracle accurate on these bodies.
  for (Oracle o : {Oracle::CG_Ball, Oracle::CG_CDHR, Oracle::CB_Bill})
  {
    checkAccurate(userBox(), o, true, 5);
    checkAccurate(longBox(5, 1000), o, true, 5);
    checkAccurate(rotatedLongBox5(), o, true, 5);
  }
  checkAccurate(longBox(10, 1000), Oracle::CB_Bill, true, 5);
  checkAccurate(longBox(10, 1000), Oracle::CG_CDHR, true, 3);
}

TEST_CASE("alt oracles: well-rounded bodies, every cooling oracle")
{
  Body c10 = cube(10);
  Body s5 = simplex(5);
  for (Oracle o : {Oracle::CG_Ball, Oracle::CG_CDHR, Oracle::CB_Bill, Oracle::CB_ABill,
                   Oracle::CB_CDHR})
  {
    checkAccurate(c10, o, false, 3);
    checkAccurate(s5, o, false, 3);
  }
}

TEST_CASE("alt oracles: CB-AcceleratedBilliard never overshoots by > 100x")
{
  // Without rounding volume_cooling_balls<AcceleratedBilliardWalk> on box5-asp100 returns
  // up to 1.7e6 x the exact volume (box10-asp1000: 1.3e14 x).  Underestimates from slow
  // mixing are a known limitation; a 10^6 overshoot is not.
  Body B = longBox(5, 100);
  for (unsigned s = 1; s <= 5; ++s)
  {
    Run r = runOracle(Oracle::CB_ABill, B, 1000u + s, false);
    INFO("seed " << s << " est=" << r.vol << " exact=" << B.exactVolume);
    CHECK(r.vol < 100 * B.exactVolume);
  }
}

// ---------------------------------------------------------------------------- slow

TEST_CASE("alt oracles: full sweep" * doctest::test_suite("slow") * doctest::skip())
{
  const int seeds = std::getenv("ALT_SEEDS") ? std::atoi(std::getenv("ALT_SEEDS")) : 10;
  const double budget = std::getenv("ALT_BUDGET") ? std::atof(std::getenv("ALT_BUDGET")) : 60;
  for (const Body& B : allBodies())
  {
    if (!selected("ALT_BODIES", B.name)) continue;
    for (bool round : {false, true})
      for (auto& [o, name] : kOracles)
      {
        if (!selected("ALT_ORACLES", label(o, round))) continue;
        // SOB is far too slow beyond n=10 at any useful accuracy.
        if (name.rfind("SOB", 0) == 0 && B.A.cols() > 10) continue;
        Cell c = evalCell(o, B, round, seeds, kTol, budget);
        report(B, o, round, c);
      }
  }
}

TEST_CASE("alt oracles: volesti inscribed_ellipsoid_rounding works in d=2")
{
  // Spectra SymEigsSolver(nev=2, ncv=min(max(10,d/5),d)) requires nev <= d-1.
  HPoly P = toHPoly(userBox());
  auto ball = P.ComputeInnerBall();
  CHECK_NOTHROW(inscribed_ellipsoid_rounding<MT, VT, NT>(P, ball.first));
}
