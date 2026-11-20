#include "volume/polytope_sampler.hpp"

#include <chrono>
#include <list>

#include "generators/boost_random_number_generator.hpp"
// uniform_billiard_walk.hpp must precede the accelerated billiard walk, which
// uses its compute_diameter without including it.
#include "random_walks/uniform_billiard_walk.hpp"
#include "random_walks/uniform_accelerated_billiard_walk.hpp"
#include "random_walks/uniform_ball_walk.hpp"
#include "random_walks/uniform_cdhr_walk.hpp"
#include "random_walks/uniform_rdhr_walk.hpp"
#include "sampling/random_point_generators.hpp"
#include "volume/gmp_walks.hpp"
#include "volume/polytope_volume.hpp"
#include "volume/sampling_policies.hpp"
#include "volume/volesti_body.hpp"

namespace ttc
{
namespace
{

using detail::HPoly;
using detail::Point;
using SampleRNG = BoostRandomNumberGenerator<boost::mt19937, double>;

double secondsSince(std::chrono::steady_clock::time_point start)
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// The body the walk runs in: P itself, or P rounded (x = shift + T y).
struct WalkBody
{
  InnerBall ball;  // of P
  PolytopeRounding rounding;

  WalkBody(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, bool round)
      : ball(certifiedInnerBall(A, b))
  {
    if (round) rounding = roundPolytope(A, b, ball);
  }
  const InnerBall& start() const { return rounding.applied ? rounding.ball : ball; }
};

// n points of hpoly (whose inner ball must be set), restarting `walk` at
// `center` for every point.
std::vector<Eigen::VectorXd> runRestartWalk(HPoly& hpoly, const Point& center, long long n,
                                            SamplerWalk walk, unsigned walkLength, SampleRNG& rng)
{
  using BilliardW = BilliardWalk::Walk<HPoly, SampleRNG>;
  using AccelBilliardW = AcceleratedBilliardWalk::Walk<HPoly, SampleRNG>;
  using BallW = BallWalk::Walk<HPoly, SampleRNG>;
  using RdhrW = RDHRWalk::Walk<HPoly, SampleRNG>;
  using CdhrW = CDHRWalk::Walk<HPoly, SampleRNG>;
  PushBackWalkPolicy policy;
  std::list<Point> samples;
  for (long long i = 0; i < n; ++i)
  {
    Point start = center;
    switch (walk)
    {
      case SamplerWalk::AcceleratedBilliard:
        RandomPointGenerator<AccelBilliardW>::apply(hpoly, start, 1, walkLength, samples, policy, rng);
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
        RandomPointGenerator<BilliardW>::apply(hpoly, start, 1, walkLength, samples, policy, rng);
        break;
    }
  }
  std::vector<Eigen::VectorXd> out;
  out.reserve(samples.size());
  for (const Point& sample : samples) out.push_back(sample.getCoefficients());
  return out;
}

}  // namespace

std::vector<Eigen::VectorXd> sampleRestartWalk(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                                               long long n, SamplerWalk walk, unsigned walkLength,
                                               unsigned seed, double& ballSeconds, bool round)
{
  const auto ballStart = std::chrono::steady_clock::now();
  const WalkBody body(A, b, round);
  const PolytopeRounding& r = body.rounding;
  HPoly hpoly = detail::makeBody(r.applied ? r.A : A, r.applied ? r.b : b, body.start());
  ballSeconds = secondsSince(ballStart);

  SampleRNG rng(static_cast<int>(A.cols()));
  rng.set_seed(seed);
  auto out = runRestartWalk(hpoly, Point(body.start().center), n, walk, walkLength, rng);
  if (r.applied)
    for (auto& y : out) y = r.shift + r.T * y;
  return out;
}

std::vector<MpVector> sampleRestartWalkGmp(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                                           long long n, SamplerWalk walk, unsigned walkLength,
                                           unsigned seed, double& ballSeconds, bool round)
{
  const auto ballStart = std::chrono::steady_clock::now();
  const WalkBody body(A, b, round);
  const PolytopeRounding& r = body.rounding;
  const Eigen::MatrixXd& Aw = r.applied ? r.A : A;
  const Eigen::VectorXd& bw = r.applied ? r.b : b;
  const InnerBall& start = body.start();
  ballSeconds = secondsSince(ballStart);

  std::vector<MpVector> out;
  switch (walk)
  {
    case SamplerWalk::Ball:
      out = sampleGmpBallWalk(Aw, bw, start.center, start.radius, n, walkLength, seed, true);
      break;
    case SamplerWalk::RDHR:
      out = sampleGmpRDHR(Aw, bw, start.center, start.radius, n, walkLength, seed, true);
      break;
    case SamplerWalk::CDHR:
      out = sampleGmpCDHR(Aw, bw, start.center, start.radius, n, walkLength, seed, true);
      break;
    case SamplerWalk::AcceleratedBilliard:
    case SamplerWalk::Billiard:
    default:
      out = sampleGmpBilliard(Aw, bw, start.center, start.radius, n, walkLength, seed, true);
      break;
  }
  if (r.applied)
  {
    const MpVector shift = r.shift.cast<MpFloat>();
    const auto T = r.T.cast<MpFloat>().eval();
    for (auto& y : out) y = shift + T * y;
  }
  return out;
}

}  // namespace ttc
