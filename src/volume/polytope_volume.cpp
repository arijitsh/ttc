#include "volume/polytope_volume.hpp"

#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <lp_lib.h>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/boost_random_number_generator.hpp"
#include "logger.hpp"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"
#include "random_walks/gaussian_ball_walk.hpp"
#include "random_walks/gaussian_cdhr_walk.hpp"
#include "random_walks/uniform_billiard_walk.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"

namespace ttc
{
namespace
{

using Kernel = Cartesian<double>;
using Point = Kernel::Point;
using HPoly = HPolytope<Point>;
// Volume RNG matches the volesti `volume` tool: a fixed compile-time seed (3),
// so the oracle estimate is deterministic.
using VolumeRNG = BoostRandomNumberGenerator<boost::mt19937, double, 3>;
// Explicitly seeded stream (oracle retries, tests).
using SeededRNG = BoostRandomNumberGenerator<boost::mt19937, double>;

using LpPtr = std::unique_ptr<lprec, decltype(&delete_lp)>;

LpPtr makeLp(int columns)
{
  LpPtr lp(make_lp(0, columns), delete_lp);
  if (!lp) throw std::runtime_error("outer-radius LP allocation failed");
  set_verbose(lp.get(), CRITICAL);
  return lp;
}

// max sign*y_j  s.t.  A y <= s, y free, on a fresh lp_solve model.
int solvePrimalExtent(const Eigen::MatrixXd& A, const Eigen::VectorXd& s, int j,
                      double sign, bool scaled, double& value)
{
  const int dim = static_cast<int>(A.cols());
  LpPtr lp = makeLp(dim);
  if (!scaled) set_scaling(lp.get(), SCALE_NONE);
  for (int k = 1; k <= dim; ++k) set_unbounded(lp.get(), k);
  std::vector<double> row(dim + 1, 0.0);
  set_add_rowmode(lp.get(), TRUE);
  for (Eigen::Index i = 0; i < A.rows(); ++i)
  {
    for (int k = 0; k < dim; ++k) row[k + 1] = A(i, k);
    if (!add_constraint(lp.get(), row.data(), LE, s(i)))
      throw std::runtime_error("outer-radius LP constraint failed");
  }
  set_add_rowmode(lp.get(), FALSE);
  std::fill(row.begin(), row.end(), 0.0);
  row[j + 1] = sign;
  if (!set_obj_fn(lp.get(), row.data()))
    throw std::runtime_error("outer-radius LP objective failed");
  set_maxim(lp.get());
  const int status = solve(lp.get());
  value = get_objective(lp.get());
  return status;
}

// The same extent via the LP dual  min s.lambda  s.t.  A^T lambda = sign*e_j,
// lambda >= 0.  It has no free variables, and any dual-feasible lambda
// certifies sign*y_j <= s.lambda on the whole polytope.
int solveDualExtent(const Eigen::MatrixXd& A, const Eigen::VectorXd& s, int j,
                    double sign, double& value)
{
  const int dim = static_cast<int>(A.cols());
  const int rows = static_cast<int>(A.rows());
  LpPtr lp = makeLp(rows);
  std::vector<double> row(rows + 1, 0.0);
  set_add_rowmode(lp.get(), TRUE);
  for (int k = 0; k < dim; ++k)
  {
    for (int i = 0; i < rows; ++i) row[i + 1] = A(i, k);
    if (!add_constraint(lp.get(), row.data(), EQ, k == j ? sign : 0.0))
      throw std::runtime_error("outer-radius LP constraint failed");
  }
  set_add_rowmode(lp.get(), FALSE);
  for (int i = 0; i < rows; ++i) row[i + 1] = s(i);
  if (!set_obj_fn(lp.get(), row.data()))
    throw std::runtime_error("outer-radius LP objective failed");
  set_minim(lp.get());
  const int status = solve(lp.get());
  if (status != OPTIMAL) return status;
  std::vector<double> raw(rows);
  get_variables(lp.get(), raw.data());
  const Eigen::VectorXd lambda = Eigen::Map<Eigen::VectorXd>(raw.data(), rows).cwiseMax(0.0);
  Eigen::VectorXd target = Eigen::VectorXd::Zero(dim);
  target(j) = sign;
  if (!((A.transpose() * lambda - target).norm() <= 1e-7)) return NUMFAILURE;
  value = s.dot(lambda);
  return status;
}

}  // namespace

// Coordinate extrema bound P by a box, hence by a ball about center.
// Like the existing inner-ball LP, these bounds use practical double arithmetic.
// The LPs are posed in y = x - center with unit-norm rows, so y = 0 is
// strictly feasible, and every extent gets a fresh model: re-optimising one
// lp_solve model under a new objective reported spurious UNBOUNDED on
// ill-conditioned zonotopes (max y_j "unbounded" while min y_j = -229 on a
// centrally symmetric parallelotope).  Rows are already equilibrated, so the
// first solve is unscaled; a non-optimal status is retried with lp_solve's
// default scaling and then through the dual before the polytope is rejected.
double outerRadius(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                   const Eigen::VectorXd& center)
{
  const int dim = static_cast<int>(A.cols());
  Eigen::MatrixXd unitA = A;
  Eigen::VectorXd slack = b - A * center;
  for (Eigen::Index i = 0; i < A.rows(); ++i)
  {
    const double norm = A.row(i).norm();
    if (norm > 0.0)
    {
      unitA.row(i) /= norm;
      slack(i) /= norm;
    }
  }
  double radius = 0.0;
  for (int j = 0; j < dim; ++j)
  {
    double extent = 0.0;
    for (const double sign : {1.0, -1.0})
    {
      double value = 0.0;
      int primal = solvePrimalExtent(unitA, slack, j, sign, false, value);
      if (primal != OPTIMAL)
        primal = solvePrimalExtent(unitA, slack, j, sign, true, value);
      if (primal != OPTIMAL)
      {
        const int dual = solveDualExtent(unitA, slack, j, sign, value);
        if (dual == INFEASIBLE && primal == UNBOUNDED)
          throw std::runtime_error("polytope is unbounded");
        if (dual != OPTIMAL)
          throw std::runtime_error("outer-radius LP failed (lp_solve status "
                                   + std::to_string(primal) + "/"
                                   + std::to_string(dual) + ")");
      }
      extent = std::max(extent, std::abs(value));
    }
    radius = std::hypot(radius, extent);
  }
  if (!std::isfinite(radius) || radius <= 0.0)
    throw std::runtime_error("invalid polytope outer radius");
  return radius;
}

namespace
{

// Chebyshev ball  max r  s.t.  a_i.x + r <= b_i  (unit-norm rows, x free,
// r >= 0) on a fresh lp_solve model; unscaled first, then default scaling.
bool solveChebyshevBall(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                        Eigen::VectorXd& center, double& radius)
{
  const int dim = static_cast<int>(A.cols());
  for (const bool scaled : {false, true})
  {
    LpPtr lp = makeLp(dim + 1);
    if (!scaled) set_scaling(lp.get(), SCALE_NONE);
    for (int k = 1; k <= dim; ++k) set_unbounded(lp.get(), k);
    std::vector<double> row(dim + 2, 0.0);
    row[dim + 1] = 1.0;
    set_add_rowmode(lp.get(), TRUE);
    for (Eigen::Index i = 0; i < A.rows(); ++i)
    {
      for (int k = 0; k < dim; ++k) row[k + 1] = A(i, k);
      if (!add_constraint(lp.get(), row.data(), LE, b(i)))
        throw std::runtime_error("Chebyshev-ball LP constraint failed");
    }
    set_add_rowmode(lp.get(), FALSE);
    std::fill(row.begin(), row.end(), 0.0);
    row[dim + 1] = 1.0;
    if (!set_obj_fn(lp.get(), row.data()))
      throw std::runtime_error("Chebyshev-ball LP objective failed");
    set_maxim(lp.get());
    if (solve(lp.get()) != OPTIMAL) continue;
    std::vector<double> values(dim + 1);
    get_variables(lp.get(), values.data());
    center = Eigen::Map<Eigen::VectorXd>(values.data(), dim);
    radius = values[dim];
    return true;
  }
  return false;
}

// volesti's ComputeInnerBall returns the max_inscribed_ball IPM iterate without
// checking that the IPM converged (it did not on 116 of 393 syntheticLRA1131
// polytopes, and stalled at a ball 11x too small on a dim-27 simplex), and
// its lp_solve fallback trusts any OPTIMAL status.  The radius is never checked
// against the polytope, yet the oracle's walk step 4r/sqrt(a*n), the clipping
// bounds and GetPrecision all scale with it.  An overestimated radius stalls the
// Gaussian walk at the centre, the annealing schedule collapses to the single
// value a = 0, and the oracle returns (pi/0)^(n/2) = inf.  An underestimated one
// slows the walk and biases the estimate low.  So re-derive the radius from the
// slacks at each candidate centre, compare against a fresh Chebyshev LP, and
// keep the larger certified ball.  A valid volesti ball within 0.1% of the LP
// optimum is returned untouched, so well-behaved inputs give bit-identical
// results.  The result is stored in P, where volesti's walks read it.
template <typename HPoly>
std::pair<typename HPoly::PointType, double> reliableInnerBall(HPoly& P)
{
  using PointT = typename HPoly::PointType;
  constexpr double kRelTol = 1e-6;       // claimed radius vs. certified slack
  constexpr double kOptimalTol = 1e-3;   // valid but clearly suboptimal ball
  const auto ball = P.ComputeInnerBall();  // also row-normalizes P
  const Eigen::MatrixXd A = P.get_mat();
  const Eigen::VectorXd b = P.get_vec();
  const double none = -std::numeric_limits<double>::infinity();
  auto slackRadius = [&](const Eigen::VectorXd& c) {
    return c.size() == A.cols() && c.allFinite() ? (b - A * c).minCoeff() : none;
  };

  const Eigen::VectorXd volestiCenter = ball.first.getCoefficients();
  const double volestiCert = slackRadius(volestiCenter);
  Eigen::VectorXd lpCenter;
  double lpRadius = 0.0;
  double lpCert = none;
  if (solveChebyshevBall(A, b, lpCenter, lpRadius))
    lpCert = std::min(lpRadius, slackRadius(lpCenter));

  const bool volestiValid = std::isfinite(ball.second) && ball.second > 0.0
                            && ball.second <= volestiCert * (1.0 + kRelTol);
  if (volestiValid && !(ball.second < lpCert * (1.0 - kOptimalTol)))
    return ball;

  std::pair<PointT, double> fixed{PointT(static_cast<unsigned>(A.cols())), -1.0};
  if (lpCert >= volestiCert && lpCert > 0.0)
    fixed = {PointT(lpCenter), lpCert};
  else if (volestiCert > 0.0)
    fixed = {ball.first, volestiCert};
  Log(2) << "inner ball: volesti radius " << ball.second << " (certified "
         << volestiCert << "), Chebyshev LP radius " << lpCert << "; using "
         << fixed.second << std::endl;
  if (fixed.second > 0.0) P.set_InnerBall(fixed);
  return fixed;
}

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

HPoly makeBody(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const InnerBall& ball)
{
  HPoly P(static_cast<unsigned>(A.cols()), A, b);
  // The state certifiedInnerBall leaves P in: rows normalised, ball stored.
  P.normalize();
  if (ball.radius > 0.0) P.set_InnerBall({Point(ball.center), ball.radius});
  return P;
}

}  // namespace

InnerBall certifiedInnerBall(const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
{
  HPoly P(static_cast<unsigned>(A.cols()), A, b);
  const auto ball = reliableInnerBall(P);
  return {ball.first.getCoefficients(), ball.second};
}

struct RoundedPolytope::Impl
{
  HPoly body;
  RoundingResult rounding;
  double innerRadius = -1.0;
};

RoundedPolytope::RoundedPolytope(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                                 const InnerBall& ball, bool round)
    : impl_(std::make_unique<Impl>())
{
  const HPoly original = makeBody(A, b, ball);
  impl_->body = original;
  impl_->innerRadius = ball.radius;
  if (!round || ball.radius <= 0.0) return;
  impl_->rounding = roundForVolume(impl_->body, Point(ball.center));
  if (!impl_->rounding.applied) return;
  const auto roundedBall = reliableInnerBall(impl_->body);
  if (roundedBall.second > 0.0 && std::isfinite(roundedBall.second))
  {
    impl_->innerRadius = roundedBall.second;
    return;
  }
  impl_->rounding = RoundingResult{};
  impl_->rounding.failure = "no inner ball after rounding";
  impl_->body = original;
}

RoundedPolytope::~RoundedPolytope() = default;
RoundedPolytope::RoundedPolytope(RoundedPolytope&&) noexcept = default;
RoundedPolytope& RoundedPolytope::operator=(RoundedPolytope&&) noexcept = default;

bool RoundedPolytope::rounded() const { return impl_->rounding.applied; }
double RoundedPolytope::det() const { return impl_->rounding.det; }
const std::string& RoundedPolytope::roundingFailure() const { return impl_->rounding.failure; }
double RoundedPolytope::innerRadius() const { return impl_->innerRadius; }

double RoundedPolytope::volume(VolumeOracle oracle, double error, unsigned walkLength) const
{
  HPoly body = impl_->body;
  VolumeRNG rng(static_cast<int>(body.dimension()));
  return impl_->rounding.det * runVolumeOracle(body, oracle, error, walkLength, rng);
}

double RoundedPolytope::volumeWithSeed(VolumeOracle oracle, double error, unsigned walkLength,
                                       unsigned seed) const
{
  HPoly body = impl_->body;
  SeededRNG rng(static_cast<int>(body.dimension()));
  rng.set_seed(seed);
  return impl_->rounding.det * runVolumeOracle(body, oracle, error, walkLength, rng);
}

}  // namespace ttc
