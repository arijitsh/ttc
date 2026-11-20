#include "volume/polytope_volume.hpp"

#include <cmath>
#include <exception>
#include <memory>
#include <string>
#include <tuple>

#include "generators/boost_random_number_generator.hpp"
#include "preprocess/inscribed_ellipsoid_rounding.hpp"
#include "random_walks/gaussian_ball_walk.hpp"
#include "random_walks/gaussian_cdhr_walk.hpp"
#include "random_walks/uniform_billiard_walk.hpp"
#include "volume/volesti_body.hpp"
#include "volume/volume_cooling_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"

namespace ttc
{
namespace
{

using detail::HPoly;
using detail::Point;
// Fixed compile-time seed 3, like the volesti `volume` tool.
using VolumeRNG = BoostRandomNumberGenerator<boost::mt19937, double, 3>;
using SeededRNG = BoostRandomNumberGenerator<boost::mt19937, double>;

struct RoundingResult
{
  double det = 1.0;     // vol(P_original) = det * vol(P_rounded)
  bool applied = false;
  std::string failure;  // why rounding was skipped (empty if applied)
  Eigen::MatrixXd T;    // x = shift + T y maps the rounded body onto P
  Eigen::VectorXd shift;
};

// Rounds P in place so that its maximum inscribed ellipsoid becomes the unit
// ball (iterated until the ellipsoid's eigenvalue ratio is <= maxEigRatio or
// maxIterations is reached), translating `interior` (a strictly interior
// point, e.g. the certified Chebyshev centre) to the origin.  If rounding
// fails or yields a non-finite or degenerate map, P is restored and det = 1.
// The caller must recompute P's inner ball afterwards.
RoundingResult roundForVolume(HPoly& P, const Point& interior, int maxIterations = 5,
                              double maxEigRatio = 6.0)
{
  RoundingResult result;
  const HPoly original = P;
  try
  {
    auto rounding = inscribed_ellipsoid_rounding<Eigen::MatrixXd, Eigen::VectorXd, double>(
        P, interior, maxIterations, maxEigRatio);
    const double det = std::get<2>(rounding);
    if (!std::isfinite(det) || det <= 0.0)
      result.failure = "non-finite or zero determinant";
    else if (!P.get_mat().allFinite() || !P.get_vec().allFinite())
      result.failure = "non-finite rounded polytope";
    else if (!(P.get_vec().array() > 0.0).all())
      result.failure = "rounding centre is not interior";
    else
    {
      result.det = det;
      result.applied = true;
      result.T = std::get<0>(rounding);
      result.shift = std::get<1>(rounding);
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
template <typename RNG>
double runVolumeOracle(HPoly& P, VolumeOracle oracle, double error, unsigned walkLength, RNG& rng)
{
  switch (oracle)
  {
    case VolumeOracle::CoolingGaussiansCdhr:
      return volume_cooling_gaussians<GaussianCDHRWalk>(P, rng, error, walkLength);
    case VolumeOracle::CoolingBallsBilliard:
      // (log vol, vol), and (-1, 0) on failure: a non-positive volume.
      return volume_cooling_balls<BilliardWalk>(P, rng, error, walkLength).second;
    case VolumeOracle::CoolingGaussiansBall:
    default:
      return volume_cooling_gaussians<GaussianBallWalk>(P, rng, error, walkLength);
  }
}

}  // namespace

PolytopeRounding roundPolytope(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                               const InnerBall& ball)
{
  PolytopeRounding out;
  if (ball.radius <= 0.0)
  {
    out.failure = "no inner ball";
    return out;
  }
  HPoly body = detail::makeBody(A, b, ball);
  const RoundingResult rounding = roundForVolume(body, Point(ball.center));
  if (!rounding.applied)
  {
    out.failure = rounding.failure;
    return out;
  }
  const auto roundedBall = detail::reliableInnerBall(body);
  if (!(roundedBall.second > 0.0 && std::isfinite(roundedBall.second)))
  {
    out.failure = "no inner ball after rounding";
    return out;
  }
  out.applied = true;
  out.det = rounding.det;
  out.T = rounding.T;
  out.shift = rounding.shift;
  out.A = body.get_mat();
  out.b = body.get_vec();
  out.ball = {roundedBall.first.getCoefficients(), roundedBall.second};
  return out;
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
  const HPoly original = detail::makeBody(A, b, ball);
  impl_->body = original;
  impl_->innerRadius = ball.radius;
  if (!round || ball.radius <= 0.0) return;
  impl_->rounding = roundForVolume(impl_->body, Point(ball.center));
  if (!impl_->rounding.applied) return;
  const auto roundedBall = detail::reliableInnerBall(impl_->body);
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
