// Rounding as the fix for ill-conditioned (long / thin / rotated) polytopes.
//
// ttc calls volume_cooling_gaussians directly on the raw H-polytope (no
// rounding).  Gaussian cooling starts from a Gaussian concentrated in the
// inner ball and walks with GaussianBallWalk, whose step radius is tied to the
// inner radius (4r/sqrt(n)); on a body with aspect ratio >> 1 it cannot mix
// within walk length 10+n^2 and the estimate collapses (box [0,2]x[0,1000]:
// 7.78 instead of 2000).
//
// volesti ships rounding preprocessors (src/volesti/include/preprocess/).  Each
// one transforms P in place into P' = T^{-1}(P - shift) and returns
// (T, shift, round_value) with
//        vol(P) = round_value * vol(P')          (round_value = |det T|)
// (see inscribed_ellipsoid_rounding.hpp: round_val *= det(L^T) with
// P.linear_transformIt(L), i.e. A <- A L, P' = L^{-1}P; svd_rounding returns
// |det T|; min_sampling_covering_ellipsoid_rounding multiplies det(L^{-1}) for
// A <- A L^{-T}.  volesti's own test/rounding_test.cpp uses
// vol = get<2>(res) * volume(P_rounded).)
//
// roundedVolume(Body, method) = round, recompute the inner ball of the rounded
// body, then run ttc's default oracle (ttcDefaultOracle, cg-cdhr).
//
// Tolerance: ttc's oracle target is epsilon/12 = 0.0667; assertions use
// 3 x 0.0667 = 0.2 relative error.
#include "doctest.h"
#include "common.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "sampling/random_point_generators.hpp"
#include "volume/sampling_policies.hpp"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"
#include "preprocess/min_sampling_covering_ellipsoid_rounding.hpp"
#include "preprocess/svd_rounding.hpp"

using namespace vt;

// ttc's default per-polytope oracle (--vol-oracle cg-cdhr): Gaussian cooling
// with coordinate hit-and-run, walk length 10 + n/10, deterministic stream.
// (common.hpp's ttcOracle is the Gaussian-ball-walk call, which after the safe
// annealing step is too slow for this suite's many oracle runs.)
static double ttcDefaultOracle(HPoly P, double epsilon = 0.8)
{
  const int n = static_cast<int>(P.dimension());
  return volume_cooling_gaussians<GaussianCDHRWalk, TtcVolumeRNG>(
      P, epsilon / 12.0, ttc::defaultOracleWalkLength(ttc::VolumeOracle::CoolingGaussiansCdhr, n));
}

namespace
{

constexpr double kOracleTarget = 0.8 / 12.0;  // ttc default epsilon / 12
constexpr double kTol = 3.0 * kOracleTarget;  // 0.2

enum class RM
{
  None,        // ttc today: no rounding
  MaxEllE,     // inscribedRoundingEigen<MAX_ELLIPSOID> (local, Spectra-free)
  LogBarrierE, // inscribedRoundingEigen<LOG_BARRIER>   (local, Spectra-free)
  MaxEll,      // inscribed_ellipsoid_rounding<MAX_ELLIPSOID>  (John ellipsoid IPM)
  LogBarrier,  // inscribed_ellipsoid_rounding<LOG_BARRIER>    (analytic center)
  VolBarrier,  // inscribed_ellipsoid_rounding<VOLUMETRIC_BARRIER>
  Vaidya,      // inscribed_ellipsoid_rounding<VAIDYA_BARRIER>
  MinCover,    // min_sampling_covering_ellipsoid_rounding<BilliardWalk> (Khachiyan)
  SVD          // svd_rounding<BilliardWalk>
};

const char* rmName(RM m)
{
  switch (m)
  {
    case RM::None: return "none";
    case RM::MaxEllE: return "maxEllE";
    case RM::LogBarrierE: return "logBarE";
    case RM::MaxEll: return "maxEll";
    case RM::LogBarrier: return "logBarrier";
    case RM::VolBarrier: return "volBarrier";
    case RM::Vaidya: return "vaidya";
    case RM::MinCover: return "minCover";
    case RM::SVD: return "svd";
  }
  return "?";
}

double now()
{
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// Axis ratio (longest / shortest semi-axis) of the maximum inscribed (John)
// ellipsoid of P, computed from P's inner-ball centre.  ~1 means well rounded.
// Used only as a diagnostic.
double johnAxisRatio(HPoly P)
{
  auto ball = P.ComputeInnerBall();
  if (!(ball.second > 0)) return NAN;
  JohnEllipsoidParams<NT> jp;
  auto [E, c, conv] = max_inscribed_ellipsoid<MT>(P.get_mat(), P.get_vec(),
                                                  ball.first.getCoefficients(), jp);
  (void)c;
  (void)conv;
  MT S = (E + E.transpose()) / 2.0;
  Eigen::SelfAdjointEigenSolver<MT> es(S);
  if (es.info() != Eigen::Success) return NAN;
  double lo = es.eigenvalues().minCoeff(), hi = es.eigenvalues().maxCoeff();
  if (!(lo > 0)) return INFINITY;
  return std::sqrt(hi / lo);
}

// Local copy of volesti's inscribed_ellipsoid_rounding loop
// (preprocess/inscribed_ellipsoid_rounding.hpp:61-129) with ONE change: the
// eigenvalues of E are computed with Eigen's dense SelfAdjointEigenSolver
// instead of Spectra.  volesti's version builds
//   Spectra::SymEigsSolver(&op, nev=2, ncv=min(max(10,d/5),d))
// which throws "nev must satisfy 1 <= nev <= n - 1" for d = 2 (line 93-94), so
// it cannot round ANY 2-D polytope (e.g. the user's box).  d is small in ttc,
// so a dense eigensolve is cheap.  Also returns the final squared axis ratio.
template <int EllType>
std::tuple<MT, VT, NT> inscribedRoundingEigen(HPoly& P, const Point& inner, int maxIter,
                                              NT maxEigRatio = 6.0, NT* lastRatio = nullptr)
{
  const int d = static_cast<int>(P.dimension());
  VT x0 = inner.getCoefficients(), center, shift = VT::Zero(d);
  MT E, L, T = MT::Identity(d, d);
  bool converged = false;
  NT round_val = 1.0, R = 100, r = 1;
  EllipsoidParams<NT> params;
  for (int iter = 1;; ++iter)
  {
    std::tie(E, center, converged) =
        compute_inscribed_ellipsoid<MT, EllType>(P.get_mat(), P.get_vec(), x0, params);
    E = (E + E.transpose()) / 2.0;
    E += MT::Identity(d, d) * 1e-8;
    Eigen::LLT<MT> llt(E.llt().solve(MT::Identity(d, d)));
    L = llt.matrixL();
    Eigen::SelfAdjointEigenSolver<MT> es(E);
    R = 1.0 / es.eigenvalues().minCoeff();
    r = 1.0 / es.eigenvalues().maxCoeff();
    P.shift(center);
    shift.noalias() += T * center;
    T.applyOnTheRight(L);
    round_val *= L.transpose().determinant();
    P.linear_transformIt(L);
    params.john_params.reg = std::max(params.john_params.reg / 10.0, 1e-10);
    P.normalize();
    x0 = VT::Zero(d);
    if ((std::abs(R / r) <= maxEigRatio && converged) || iter >= maxIter) break;
  }
  if (lastRatio) *lastRatio = R / r;
  return {T, shift, std::abs(round_val)};
}

struct Rounded
{
  HPoly P;
  double roundVal = 1.0;
  double secs = 0.0;
  bool ok = true;
  double eigRatio = NAN;  // last (R/r) seen by the Eigen rounding loop
  std::string err;
};

// Round B with method m.  maxIter = rounding iterations for the inscribed
// ellipsoid routines (they stop early once eig ratio <= 6 and converged).
Rounded roundBody(const Body& B, RM m, int maxIter = 20, unsigned seed = 7)
{
  Rounded R{toHPoly(B)};
  HPoly& P = R.P;
  auto ball = P.ComputeInnerBall();  // normalizes P; same role as ttc's reliableInnerBall
  if (!(ball.second > 0))
  {
    R.ok = false;
    return R;
  }
  const int n = static_cast<int>(P.dimension());
  const double t0 = now();
  std::tuple<MT, VT, NT> res;
  TtcSampleRNG rng(n);
  rng.set_seed(seed);
  try {
  switch (m)
  {
    case RM::None: res = {MT::Identity(n, n), VT::Zero(n), 1.0}; break;
    case RM::MaxEllE:
      res = inscribedRoundingEigen<EllipsoidType::MAX_ELLIPSOID>(P, ball.first, maxIter, 6.0,
                                                                 &R.eigRatio);
      break;
    case RM::LogBarrierE:
      res = inscribedRoundingEigen<EllipsoidType::LOG_BARRIER>(P, ball.first, maxIter, 6.0,
                                                               &R.eigRatio);
      break;
    case RM::MaxEll:
      res = inscribed_ellipsoid_rounding<MT, VT, NT, HPoly, Point, EllipsoidType::MAX_ELLIPSOID>(
          P, ball.first, maxIter);
      break;
    case RM::LogBarrier:
      res = inscribed_ellipsoid_rounding<MT, VT, NT, HPoly, Point, EllipsoidType::LOG_BARRIER>(
          P, ball.first, maxIter);
      break;
    case RM::VolBarrier:
      res = inscribed_ellipsoid_rounding<MT, VT, NT, HPoly, Point,
                                         EllipsoidType::VOLUMETRIC_BARRIER>(P, ball.first,
                                                                            maxIter);
      break;
    case RM::Vaidya:
      res = inscribed_ellipsoid_rounding<MT, VT, NT, HPoly, Point, EllipsoidType::VAIDYA_BARRIER>(
          P, ball.first, maxIter);
      break;
    case RM::MinCover:
      res = min_sampling_covering_ellipsoid_rounding<BilliardWalk, MT, VT>(P, ball, 10 + 10 * n,
                                                                             rng);
      break;
    case RM::SVD: res = svd_rounding<BilliardWalk, MT, VT>(P, ball, 10 + 10 * n, rng); break;
  }
  } catch (const std::exception& e) {
    R.secs = now() - t0;
    R.ok = false;
    R.err = e.what();
    return R;
  }
  R.secs = now() - t0;
  R.roundVal = std::get<2>(res);
  if (!std::isfinite(R.roundVal) || R.roundVal <= 0) R.ok = false;
  // Recompute the inner ball of the rounded body (transforms reset has_ball).
  auto b2 = P.ComputeInnerBall();
  if (!(b2.second > 0) || !std::isfinite(b2.second)) R.ok = false;
  return R;
}

struct VolRun
{
  double vol = NAN, roundSecs = 0, oracleSecs = 0, roundVal = 1;
  bool ok = false;
  double rel(const Body& B) const { return vol / B.exactVolume - 1.0; }
};

// round + ttc-configured oracle.
VolRun roundedVolume(const Body& B, RM m = RM::MaxEll, int maxIter = 20)
{
  VolRun V;
  Rounded R = roundBody(B, m, maxIter);
  V.roundSecs = R.secs;
  V.roundVal = R.roundVal;
  if (!R.ok)
  {
    if (!R.err.empty()) std::printf("  [%s threw: %s]\n", rmName(m), R.err.c_str());
    return V;
  }
  const double t0 = now();
  V.vol = R.roundVal * ttcDefaultOracle(R.P);
  V.oracleSecs = now() - t0;
  V.ok = std::isfinite(V.vol) && V.vol > 0;
  return V;
}

void report(const Body& B, RM m, const VolRun& V)
{
  std::printf("  %-34s %-10s est=%-12.6g exact=%-12.6g rel=%+8.4f round=%.3fs oracle=%.3fs\n",
              B.name.c_str(), rmName(m), V.vol, B.exactVolume, V.rel(B), V.roundSecs,
              V.oracleSecs);
  std::fflush(stdout);
}

// ------------------------------------------------------------ test bodies
Body longBox(int n, double aspect)
{
  std::vector<double> lo(n, 0.0), hi(n, 1.0);
  hi[n - 1] = aspect;
  Body B = box(lo, hi);
  char buf[64];
  std::snprintf(buf, sizeof buf, "longbox%d[a=%g]", n, aspect);
  B.name = buf;
  return B;
}

Body rotatedLongBox(int n, double aspect, unsigned seed)
{
  Body B = longBox(n, aspect);
  VT t = VT::Constant(n, 3.0);
  return linearImage(B, randomRotation(n, seed), t, "rot" + std::to_string(seed));
}

// Standard simplex stretched by `aspect` along the last axis, then rotated.
Body longSimplex(int n, double aspect, bool rotate)
{
  Body S = simplex(n);
  MT D = MT::Identity(n, n);
  D(n - 1, n - 1) = aspect;
  Body L = linearImage(S, D, VT::Zero(n), "x" + std::to_string(int(aspect)));
  if (rotate) L = linearImage(L, randomRotation(n, 11), VT::Constant(n, -2.0), "rot");
  return L;
}

Body slab()
{
  Body B = box({0, 0, 0, 0, 0}, {1e-4, 1, 1, 1, 1});
  B.name = "slab5[1e-4]";
  return B;
}

std::vector<Body> hardSuite()
{
  std::vector<Body> v;
  Body lb = box({0, 0}, {2, 1000});
  lb.name = "box[0,2]x[0,1000]";
  v.push_back(lb);
  for (int n : {2, 5, 10})
    for (double a : {10.0, 100.0, 1e3, 1e4}) v.push_back(longBox(n, a));
  for (int n : {2, 5, 10}) v.push_back(rotatedLongBox(n, 1e3, 5));
  for (int n : {2, 5, 10}) v.push_back(longSimplex(n, 1e3, n > 2));
  v.push_back(slab());
  return v;
}

}  // namespace

// Method used by the asserting tests: volesti's public API, exactly the call ttc
// should make:
//   inscribed_ellipsoid_rounding<MT, VT, NT, HPoly, Point, EllipsoidType::MAX_ELLIPSOID>
//       (P, innerBall.first, /*max_iterations=*/20, /*max_eig_ratio=*/6.0)
// then vol = round_val * volume_cooling_gaussians(P_rounded, eps/12, 10+n^2).
// NOTE: before commit 6494cf2 this threw from Spectra in d = 2 ("nev must
// satisfy 1 <= nev <= n - 1") and, for an already-round body such as [0,1]^10,
// "TridiagEigen: eigen decomposition failed" (max_inscribed_ellipsoid.hpp ~131,
// inscribed_ellipsoid_rounding.hpp ~93).  Both are guarded by the dense Eigen
// fallback now; the d = 1,2 test and longbox10[a=1] / cube cases check that.
constexpr RM kRec = RM::MaxEll;

TEST_CASE("rounding: volesti inscribed_ellipsoid_rounding must not throw in d = 1, 2")
{
  // Spectra::SymEigsSolver(nev = 2) requires nev <= n - 1.
  for (int n : {1, 2, 3})
    for (RM m : {RM::MaxEll, RM::LogBarrier})
    {
      Body B = longBox(n, 50.0);
      Rounded R = roundBody(B, m);
      CAPTURE(n);
      CAPTURE(rmName(m));
      CAPTURE(R.err);
      CHECK(R.err.empty());
      CHECK(R.ok);
    }
}

// =====================================================================
// Recommended pipeline: max-inscribed-ellipsoid rounding + ttc oracle.
// =====================================================================

TEST_CASE("rounding: user-reported long box [0,2]x[0,1000]")
{
  Body B = box({0, 0}, {2, 1000});
  B.name = "box[0,2]x[0,1000]";
  VolRun raw = roundedVolume(B, RM::None);
  VolRun rnd = roundedVolume(B, kRec);
  report(B, RM::None, raw);
  report(B, kRec, rnd);
  CHECK(rnd.ok);
  CHECK(std::abs(rnd.rel(B)) < kTol);
}

TEST_CASE("rounding: axis-aligned aspect sweep 1..1e4, dims 2,5,10")
{
  for (int n : {2, 5, 10})
    for (double a : {1.0, 10.0, 100.0, 1e3, 1e4})
    {
      Body B = longBox(n, a);
      VolRun V = roundedVolume(B, kRec);
      report(B, kRec, V);
      CAPTURE(B.name);
      CHECK(V.ok);
      CHECK(std::abs(V.rel(B)) < kTol);
    }
}

TEST_CASE("rounding: rotated + shifted long boxes (aspect 1e3)")
{
  for (int n : {2, 5, 10})
    for (unsigned s : {5u, 6u})
    {
      Body B = rotatedLongBox(n, 1e3, s);
      VolRun V = roundedVolume(B, kRec);
      report(B, kRec, V);
      CAPTURE(B.name);
      CHECK(V.ok);
      CHECK(std::abs(V.rel(B)) < kTol);
    }
}

TEST_CASE("rounding: long simplices (stretched 1e3, rotated)")
{
  for (int n : {2, 5, 10})
  {
    Body B = longSimplex(n, 1e3, n > 2);
    VolRun V = roundedVolume(B, kRec);
    report(B, kRec, V);
    CAPTURE(B.name);
    CHECK(V.ok);
    CHECK(std::abs(V.rel(B)) < kTol);
  }
}

TEST_CASE("rounding: near-degenerate slab [0,1e-4]x[0,1]^4")
{
  Body B = slab();
  VolRun V = roundedVolume(B, kRec);
  report(B, kRec, V);
  CHECK(V.ok);
  CHECK(std::abs(V.rel(B)) < kTol);
}

// The unrounded comparison is informational only, and with the safe
// annealing step the raw CG-ball oracle on this aspect-1e4 slab takes far
// longer than the default suite's budget (killed after 1500 s).
TEST_CASE("rounding: near-degenerate slab, unrounded comparison" * doctest::test_suite("slow"))
{
  Body B = slab();
  VolRun raw = roundedVolume(B, RM::None);
  report(B, RM::None, raw);
}

TEST_CASE("rounding: does not hurt well-rounded cubes n=2..10")
{
  for (int n = 2; n <= 10; ++n)
  {
    Body B = cube(n);
    VolRun V = roundedVolume(B, kRec);
    report(B, kRec, V);
    CAPTURE(B.name);
    CHECK(V.ok);
    CHECK(std::abs(V.rel(B)) < kTol);
  }
}

TEST_CASE("rounding: round value convention vol(P) = round_val * vol(P')"
          " checked exactly on a linear image of a cube")
{
  // P = Q [-1,1]^n with Q ill-conditioned.  After rounding, P' must have
  // volume vol(P)/round_val; compare with a high-accuracy oracle on P'.
  const int n = 4;
  MT Q = randomRotation(n, 3);
  Q.col(0) *= 200.0;
  Body B = linearImage(cube(n), Q, VT::Zero(n), "skew");
  for (RM m : {kRec, RM::LogBarrierE})
  {
    Rounded R = roundBody(B, m);
    REQUIRE(R.ok);
    std::vector<double> est;
    for (unsigned s = 1; s <= 4; ++s)
      est.push_back(R.roundVal * gaussianOracle(R.P, s, 0.03, 10 + n * n));
    Stats S = relStats(est, B.exactVolume, kTol);
    std::printf("  %s %s roundVal=%g %s\n", B.name.c_str(), rmName(m), R.roundVal,
                fmt(S).c_str());
    CHECK(S.nonFinite == 0);
    CHECK(std::abs(S.meanRel) < 0.1);
  }
}

// Statistical check of the recommended pipeline across oracle seeds.
TEST_CASE("rounding: long box across 8 oracle seeds")
{
  Body B = rotatedLongBox(5, 1e3, 9);
  Rounded R = roundBody(B, kRec);
  REQUIRE(R.ok);
  std::vector<double> est;
  for (unsigned s = 1; s <= 8; ++s)
    est.push_back(R.roundVal * gaussianOracle(R.P, s, kOracleTarget, ttcWalkLength(5)));
  Stats S = relStats(est, B.exactVolume, kTol);
  std::printf("  %s maxEll x8 seeds: %s\n", B.name.c_str(), fmt(S).c_str());
  CHECK(S.nonFinite == 0);
  CHECK(S.fracWithin >= 0.875);
  CHECK(std::abs(S.meanRel) < 2 * kOracleTarget);
}

// =====================================================================
// Slow comparisons (run with -ts=slow or --no-skip).
// =====================================================================

TEST_CASE("rounding-compare: every method on the hard suite" * doctest::test_suite("slow")
          * doctest::skip())
{
  // Optional filters: ROUND_BODY=<substring of body name>, ROUND_METHOD=<rmName>.
  const char* fb = std::getenv("ROUND_BODY");
  const char* fm = std::getenv("ROUND_METHOD");
  for (const Body& B : hardSuite())
  {
    if (fb && B.name.find(fb) == std::string::npos) continue;
    for (RM m : {RM::None, RM::MaxEllE, RM::LogBarrierE, RM::MaxEll, RM::LogBarrier,
                 RM::VolBarrier, RM::Vaidya, RM::MinCover, RM::SVD})
    {
      if (fm && std::string(fm) != rmName(m)) continue;
      VolRun V = roundedVolume(B, m);
      report(B, m, V);
    }
    std::printf("\n");
  }
}

TEST_CASE("rounding-compare: John-ellipsoid axis ratio vs rounding iterations"
          * doctest::test_suite("slow") * doctest::skip())
{
  for (const Body& B : {longBox(10, 1e4), rotatedLongBox(10, 1e3, 5), longSimplex(10, 1e3, true),
                        slab()})
    for (RM m : {RM::MaxEllE, RM::LogBarrierE})
      for (int it : {1, 2, 3, 5, 20})
      {
        Rounded R = roundBody(B, m, it);
        double ar = R.ok ? johnAxisRatio(R.P) : NAN;
        std::printf("  %-34s %-10s iters<=%-2d roundVal=%-12.6g axisRatio=%-10.4g eigRatio=%-10.4g %.3fs\n",
                    B.name.c_str(), rmName(m), it, R.roundVal, ar, R.eigRatio, R.secs);
        std::fflush(stdout);
      }
}

TEST_CASE("rounding-compare: cubes n=2..15, all methods" * doctest::test_suite("slow")
          * doctest::skip())
{
  for (int n = 2; n <= 15; ++n)
  {
    Body B = cube(n);
    for (RM m : {RM::None, RM::MaxEllE, RM::LogBarrierE, RM::MinCover, RM::SVD})
      report(B, m, roundedVolume(B, m));
  }
}
