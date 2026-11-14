// Shared helpers for the standalone volesti unit tests (tests/volesti/).
//
// Each test_*.cpp in this directory is compiled into its own doctest binary by
// tests/volesti/Makefile.  The helpers below build H-polytopes with known
// exact volumes and run the volesti volume oracle in exactly the configuration
// ttc's LRA engine uses (src/volume/lra_volume.cpp):
//
//   volume_cooling_gaussians<GaussianBallWalk,
//                            BoostRandomNumberGenerator<boost::mt19937,double,3>>
//       (hpoly, epsilon / 12, 10 + dim^2)
//
// with the inner ball pre-set on the polytope.
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <boost/random.hpp>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/boost_random_number_generator.hpp"
#include "random_walks/gaussian_ball_walk.hpp"
#include "random_walks/gaussian_cdhr_walk.hpp"
#include "random_walks/uniform_billiard_walk.hpp"
#include "random_walks/uniform_ball_walk.hpp"
#include "random_walks/uniform_cdhr_walk.hpp"
#include "random_walks/uniform_rdhr_walk.hpp"
#include "random_walks/uniform_accelerated_billiard_walk.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/polytope_sampler.hpp"  // src/: ttc's sampler
#include "volume/polytope_volume.hpp"   // src/: ttc's inner ball, rounding + oracle

namespace vt
{

using NT = double;
using Kernel = Cartesian<NT>;
using Point = Kernel::Point;
using HPoly = HPolytope<Point>;
using MT = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>;
using VT = Eigen::Matrix<NT, Eigen::Dynamic, 1>;

// RNGs exactly as ttc declares them.
using TtcVolumeRNG = BoostRandomNumberGenerator<boost::mt19937, double, 3>;  // fixed seed 3
using TtcSampleRNG = BoostRandomNumberGenerator<boost::mt19937, double>;     // clock seeded

// ---------------------------------------------------------------- polytopes
struct Body
{
  std::string name;
  MT A;
  VT b;
  double exactVolume = 0.0;  // 0 if unknown
};

inline HPoly toHPoly(const Body& B)
{
  return HPoly(static_cast<unsigned>(B.A.cols()), B.A, B.b);
}

// Axis-aligned box prod [lo_i, hi_i].
inline Body box(const std::vector<double>& lo, const std::vector<double>& hi)
{
  const int n = static_cast<int>(lo.size());
  Body B;
  B.name = "box" + std::to_string(n);
  B.A = MT::Zero(2 * n, n);
  B.b = VT::Zero(2 * n);
  double vol = 1.0;
  for (int i = 0; i < n; ++i)
  {
    B.A(2 * i, i) = 1.0;
    B.b(2 * i) = hi[i];
    B.A(2 * i + 1, i) = -1.0;
    B.b(2 * i + 1) = -lo[i];
    vol *= hi[i] - lo[i];
  }
  B.exactVolume = vol;
  return B;
}

// Cube [lo, hi]^n.
inline Body cube(int n, double lo = -1.0, double hi = 1.0)
{
  Body B = box(std::vector<double>(n, lo), std::vector<double>(n, hi));
  B.name = "cube" + std::to_string(n) + "[" + std::to_string(lo) + "," + std::to_string(hi) + "]";
  return B;
}

// Standard simplex {x >= 0, sum x <= s}: volume s^n / n!.
inline Body simplex(int n, double s = 1.0)
{
  Body B;
  B.name = "simplex" + std::to_string(n);
  B.A = MT::Zero(n + 1, n);
  B.b = VT::Zero(n + 1);
  for (int i = 0; i < n; ++i) B.A(i, i) = -1.0;
  B.A.row(n).setOnes();
  B.b(n) = s;
  B.exactVolume = std::pow(s, n) / std::tgamma(n + 1.0);
  return B;
}

// Cross-polytope {sum |x_i| <= 1}: 2^n facets, volume 2^n / n!.
inline Body crossPolytope(int n)
{
  Body B;
  B.name = "cross" + std::to_string(n);
  const int m = 1 << n;
  B.A = MT::Zero(m, n);
  B.b = VT::Ones(m);
  for (int k = 0; k < m; ++k)
    for (int i = 0; i < n; ++i) B.A(k, i) = (k >> i) & 1 ? 1.0 : -1.0;
  B.exactVolume = std::pow(2.0, n) / std::tgamma(n + 1.0);
  return B;
}

// x -> Q x + t applied to body (volume scales by |det Q|).
inline Body linearImage(const Body& B, const MT& Q, const VT& t, const std::string& tag)
{
  // y = Qx + t  <=>  x = Q^{-1}(y - t);  A x <= b  <=>  A Q^{-1} y <= b + A Q^{-1} t
  Body R;
  R.name = B.name + "+" + tag;
  MT Qinv = Q.inverse();
  R.A = B.A * Qinv;
  R.b = B.b + R.A * t;
  R.exactVolume = B.exactVolume * std::abs(Q.determinant());
  return R;
}

inline MT randomRotation(int n, unsigned seed)
{
  std::mt19937 g(seed);
  std::normal_distribution<double> N(0, 1);
  MT M(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) M(i, j) = N(g);
  Eigen::HouseholderQR<MT> qr(M);
  MT Q = qr.householderQ();
  return Q;
}

// --------------------------------------------------------------- inner ball
// Set the exact Chebyshev ball (center, radius) on P (after normalizing rows),
// mirroring what ttc does with reliableInnerBall.
inline void setBall(HPoly& P, const VT& center, double radius)
{
  P.normalize();
  P.set_InnerBall({Point(center), radius});
}

// Let volesti compute its own inner ball (normalizes P).
inline std::pair<Point, double> volestiBall(HPoly& P) { return P.ComputeInnerBall(); }

// ------------------------------------------------------------------ oracles
inline unsigned ttcWalkLength(int n) { return 10u + static_cast<unsigned>(n * n); }

// The raw volesti oracle exactly as ttc calls it, on whatever body it is given
// (deterministic: VolumeRNG has fixed seed 3).  No rounding: ttc rounds first,
// see ttcVolume().
inline double ttcOracle(HPoly P, double epsilon = 0.8)
{
  const int n = static_cast<int>(P.dimension());
  return volume_cooling_gaussians<GaussianBallWalk, TtcVolumeRNG>(P, epsilon / 12.0,
                                                                  ttcWalkLength(n));
}

// ttc's full per-polytope pipeline (src/volume/lra_volume.cpp phase 1b),
// exactly as ttc runs it: certified inner ball, max-inscribed-ellipsoid
// rounding, the selected oracle, times |det T|.
inline ttc::RoundedPolytope ttcBody(const MT& A, const VT& b, bool round = true)
{
  return ttc::RoundedPolytope(A, b, ttc::certifiedInnerBall(A, b), round);
}
inline double ttcVolume(const HPoly& P, double epsilon = 0.8,
                        ttc::VolumeOracle oracle = ttc::kDefaultVolumeOracle)
{
  const int n = static_cast<int>(P.dimension());
  return ttcBody(P.get_mat(), P.get_vec())
      .volume(oracle, epsilon / 12.0, ttc::defaultOracleWalkLength(oracle, n));
}
// ttcVolume on an explicitly seeded stream (for statistics across seeds).
inline double ttcVolumeSeeded(const HPoly& P, unsigned seed, double epsilon = 0.8,
                              ttc::VolumeOracle oracle = ttc::kDefaultVolumeOracle)
{
  const int n = static_cast<int>(P.dimension());
  return ttcBody(P.get_mat(), P.get_vec())
      .volumeWithSeed(oracle, epsilon / 12.0, ttc::defaultOracleWalkLength(oracle, n), seed);
}

// ------------------------------------------------------------------ sampling
// ttc's single-polytope sampler (src/volume/polytope_sampler.cpp), exactly as
// lra_volume.cpp calls it.
//   ttcProposals: n points, each the endpoint of a fresh walk of walkLength
//                 steps started at the certified Chebyshev centre (default
//                 walk length 10 + n^2, ttc's lraDefaultWalkLength).
//   ttcLatticeSamples: GenerateSamples -- proposals from the facet-inflated
//                 body, stochastic rounding to 10^-b Z^n, rejection against P.
inline std::vector<VT> ttcProposals(const MT& A, const VT& b, long long n, unsigned seed,
                                    ttc::SamplerWalk walk = ttc::SamplerWalk::Billiard,
                                    int walkLength = -1)
{
  const unsigned wl =
      walkLength >= 0 ? unsigned(walkLength) : ttc::defaultSampleWalkLength(walk, A.cols());
  double ballSeconds = 0;
  return ttc::sampleRestartWalk(A, b, n, walk, wl, seed, ballSeconds);
}

inline std::vector<VT> ttcLatticeSamples(const MT& A, const VT& b, long long n, int precision,
                                         unsigned seed,
                                         ttc::SamplerWalk walk = ttc::SamplerWalk::Billiard,
                                         double epsilon = 0.8, double delta = 0.2,
                                         int walkLength = -1)
{
  const unsigned wl =
      walkLength >= 0 ? unsigned(walkLength) : ttc::defaultSampleWalkLength(walk, A.cols());
  unsigned walkSeed = seed * 2654435761u + 1u;
  auto generate = [&](const MT& Ai, const VT& bi, long long k, double& t) {
    double ballSeconds = 0;
    t = 0;
    return ttc::sampleRestartWalk(Ai, bi, k, walk, wl, walkSeed++, ballSeconds);
  };
  double samplingTime = 0;
  return ttc::sampleLattice<double>(A, b, n, precision, epsilon / 12.0, delta, seed, generate,
                                    samplingTime);
}

// Same oracle on an explicitly seeded RNG (for statistics across seeds).
inline double gaussianOracle(HPoly P, unsigned seed, double error, unsigned walkLength)
{
  TtcSampleRNG rng(static_cast<int>(P.dimension()));
  rng.set_seed(seed);
  return volume_cooling_gaussians<GaussianBallWalk>(P, rng, error, walkLength);
}

// ---------------------------------------------------------------- statistics
struct Stats
{
  double mean = 0, sd = 0, minRel = 0, maxRel = 0, meanRel = 0;
  int n = 0, nonFinite = 0;
  double fracWithin = 0;  // fraction with |rel err| <= tol
};

// rel error = estimate/exact - 1.
inline Stats relStats(const std::vector<double>& est, double exact, double tol)
{
  Stats S;
  S.n = static_cast<int>(est.size());
  std::vector<double> rel;
  for (double e : est)
  {
    if (!std::isfinite(e) || e <= 0) { ++S.nonFinite; continue; }
    rel.push_back(e / exact - 1.0);
  }
  if (rel.empty()) return S;
  S.meanRel = std::accumulate(rel.begin(), rel.end(), 0.0) / rel.size();
  S.mean = exact * (1.0 + S.meanRel);
  double v = 0;
  for (double r : rel) v += (r - S.meanRel) * (r - S.meanRel);
  S.sd = rel.size() > 1 ? std::sqrt(v / (rel.size() - 1)) : 0.0;
  S.minRel = *std::min_element(rel.begin(), rel.end());
  S.maxRel = *std::max_element(rel.begin(), rel.end());
  S.fracWithin = double(std::count_if(rel.begin(), rel.end(),
                                      [&](double r) { return std::abs(r) <= tol; }))
                 / S.n;
  return S;
}

inline std::string fmt(const Stats& S)
{
  char buf[256];
  std::snprintf(buf, sizeof buf,
                "meanRel=%+.4f sd=%.4f min=%+.4f max=%+.4f within=%.2f nonfinite=%d/%d",
                S.meanRel, S.sd, S.minRel, S.maxRel, S.fracWithin, S.nonFinite, S.n);
  return buf;
}

}  // namespace vt
