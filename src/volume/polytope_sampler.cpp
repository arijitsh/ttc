#include "volume/polytope_sampler.hpp"

#include <chrono>
#include <list>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/boost_random_number_generator.hpp"
// uniform_billiard_walk.hpp must precede the accelerated billiard walk, which
// uses its compute_diameter without including it.
#include "random_walks/uniform_billiard_walk.hpp"
#include "random_walks/uniform_accelerated_billiard_walk.hpp"
#include "random_walks/uniform_ball_walk.hpp"
#include "random_walks/uniform_cdhr_walk.hpp"
#include "random_walks/uniform_rdhr_walk.hpp"
#include "sampling/random_point_generators.hpp"
#include "volume/polytope_volume.hpp"
#include "volume/sampling_policies.hpp"

namespace ttc
{
namespace
{

using Kernel = Cartesian<double>;
using Point = Kernel::Point;
using HPoly = HPolytope<Point>;
// Sampling RNG matches the volesti `sample` tool: explicitly seeded.
using SampleRNG = BoostRandomNumberGenerator<boost::mt19937, double>;

// Draws n points from hpoly (whose inner ball must be set), restarting the
// selected walk at `center` for every point.
template <typename HPoly, typename RNG>
std::vector<Eigen::VectorXd> runRestartWalk(HPoly& hpoly,
                                               const typename HPoly::PointType& center,
                                               long long n, SamplerWalk walk,
                                               unsigned walkLength, RNG& rng)
{
  using Point = typename HPoly::PointType;
  using BilliardW = typename BilliardWalk::template Walk<HPoly, RNG>;
  using AccelBilliardW = typename AcceleratedBilliardWalk::template Walk<HPoly, RNG>;
  using BallW = typename BallWalk::template Walk<HPoly, RNG>;
  using RdhrW = typename RDHRWalk::template Walk<HPoly, RNG>;
  using CdhrW = typename CDHRWalk::template Walk<HPoly, RNG>;
  PushBackWalkPolicy policy;
  std::list<Point> samples;
  for (long long i = 0; i < n; ++i)
  {
    Point start = center;
    switch (walk)
    {
      case SamplerWalk::AcceleratedBilliard:
        RandomPointGenerator<AccelBilliardW>::apply(hpoly, start, 1, walkLength, samples,
                                                    policy, rng);
        break;
      case SamplerWalk::Ball:
        RandomPointGenerator<BallW>::apply(hpoly, start, 1, walkLength, samples, policy, rng);
        break;
      case SamplerWalk::RDHR:
        RandomPointGenerator<RdhrW>::apply(hpoly, start, 1, walkLength, samples, policy, rng);
        break;
      case SamplerWalk::CDHR:
        RandomPointGenerator<CdhrW>::apply(hpoly, start, 1, walkLength, samples, policy, rng);
        break;
      case SamplerWalk::Billiard:
      default:
        RandomPointGenerator<BilliardW>::apply(hpoly, start, 1, walkLength, samples, policy,
                                               rng);
        break;
    }
  }
  std::vector<Eigen::VectorXd> out;
  out.reserve(samples.size());
  for (const Point& sample : samples) out.push_back(sample.getCoefficients());
  return out;
}

}  // namespace

std::vector<Eigen::VectorXd> sampleRestartWalk(const Eigen::MatrixXd& A,
                                               const Eigen::VectorXd& b, long long n,
                                               SamplerWalk walk, unsigned walkLength,
                                               unsigned seed, double& ballSeconds, bool round)
{
  const auto ballStart = std::chrono::steady_clock::now();
  const InnerBall ball = certifiedInnerBall(A, b);
  PolytopeRounding rounding;
  if (round) rounding = roundPolytope(A, b, ball);
  const Eigen::MatrixXd& Aw = rounding.applied ? rounding.A : A;
  const Eigen::VectorXd& bw = rounding.applied ? rounding.b : b;
  const InnerBall& start = rounding.applied ? rounding.ball : ball;
  HPoly hpoly(static_cast<unsigned>(A.cols()), Aw, bw);
  hpoly.normalize();
  if (start.radius > 0.0) hpoly.set_InnerBall({Point(start.center), start.radius});
  ballSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - ballStart).count();
  SampleRNG rng(static_cast<int>(A.cols()));
  rng.set_seed(seed);
  auto out = runRestartWalk(hpoly, Point(start.center), n, walk, walkLength, rng);
  if (rounding.applied)
    for (auto& y : out) y = rounding.shift + rounding.T * y;
  return out;
}

}  // namespace ttc
