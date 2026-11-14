#pragma once

// Single-polytope volume oracle used by the LRA engine (lra_volume.cpp) and by
// the standalone tests in tests/volesti/.
//
// Cousins & Vempala's practical volume algorithm (Math. Prog. Comp. 2016,
// Fig. 2) starts with `T = Round(K)` and returns |det T| * vol(T^-1 K): the
// random walks inside every volesti volume routine take steps proportional to
// the inner radius, so on a skinny body (R/r >> sqrt n) they never cross it
// within their step budget and the estimate collapses (box [0,2]x[0,1000]:
// 7.78 instead of 2000).  None of volesti's volume routines rounds internally,
// so roundForVolume() applies the max-inscribed-ellipsoid rounding of
// Cousins-Vempala's preprocess.m (volesti's inscribed_ellipsoid_rounding)
// first.  See text-for-agents/ttc-volume-vs-theory.md.

#include <cmath>
#include <exception>
#include <string>

#include <Eigen/Dense>

#include "preprocess/inscribed_ellipsoid_rounding.hpp"
#include "random_walks/gaussian_ball_walk.hpp"
#include "random_walks/gaussian_cdhr_walk.hpp"
#include "random_walks/uniform_billiard_walk.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_oracle_kind.hpp"

namespace ttc
{

struct RoundingResult
{
  double det = 1.0;     // vol(P_original) = det * vol(P_rounded)
  bool applied = false;
  std::string failure;  // why rounding was skipped (empty if applied)
};

// Rounds P in place so that its maximum inscribed ellipsoid becomes the unit
// ball (iterated until the ellipsoid's eigenvalue ratio is <= maxEigRatio or
// maxIterations is reached), translating `interior` (a strictly interior
// point, e.g. the certified Chebyshev centre) to the origin.  If rounding
// fails or yields a non-finite / degenerate map, P is restored and det = 1.
// The caller must recompute P's inner ball afterwards.
template <typename HPoly>
RoundingResult roundForVolume(HPoly& P, const typename HPoly::PointType& interior,
                              int maxIterations = 5, double maxEigRatio = 6.0)
{
  using MT = Eigen::MatrixXd;
  using VT = Eigen::VectorXd;
  RoundingResult result;
  const HPoly original = P;
  try
  {
    auto rounding = inscribed_ellipsoid_rounding<MT, VT, double>(P, interior, maxIterations,
                                                                 maxEigRatio);
    const double det = std::get<2>(rounding);
    const MT& A = P.get_mat();
    const VT& b = P.get_vec();
    if (!std::isfinite(det) || det <= 0.0)
      result.failure = "non-finite or zero determinant";
    else if (!A.allFinite() || !b.allFinite())
      result.failure = "non-finite rounded polytope";
    else if (!(b.array() > 0.0).all())
      result.failure = "rounding centre is not interior";
    else
    {
      result.det = det;
      result.applied = true;
      return result;
    }
  }
  catch (const std::exception& e)
  {
    result.failure = e.what();
  }
  P = original;
  return result;
}

// Runs the selected oracle on P, whose inner ball must already be set.
// Returns a non-finite or non-positive value on failure.
template <typename HPoly, typename RNG>
double runVolumeOracle(HPoly& P, VolumeOracle oracle, double error, unsigned walkLength,
                       RNG& rng)
{
  switch (oracle)
  {
    case VolumeOracle::CoolingGaussiansCdhr:
      return volume_cooling_gaussians<GaussianCDHRWalk>(P, rng, error, walkLength);
    case VolumeOracle::CoolingBallsBilliard:
    {
      // cooling balls returns (log vol, vol), and (-1, 0) on failure, which
      // callers see as a non-positive volume.
      const auto logAndVol = volume_cooling_balls<BilliardWalk>(P, rng, error, walkLength);
      return logAndVol.second;
    }
    case VolumeOracle::CoolingGaussiansBall:
    default:
      return volume_cooling_gaussians<GaussianBallWalk>(P, rng, error, walkLength);
  }
}

}  // namespace ttc
