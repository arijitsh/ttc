#pragma once

// Uniform sampling from one polytope, as used by the LRA engine's union
// algorithm (lra_volume.cpp) and by the standalone tests in tests/volesti/.
//
//   sampleRestartWalk: n proposals, each from a fresh run of the selected
//     volesti walk started at the Chebyshev centre (walkLength steps, endpoint
//     kept); volesti stays inside polytope_sampler.cpp.
//   sampleLattice:     GenerateSamples of the manuscript -- inflate the facets
//     by 10^-b * kappa * |A_j|, draw proposals from the inflated body, round
//     each coordinate stochastically to 10^-b Z, and reject against P.

#include <chrono>
#include <cmath>
#include <random>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>

#include "volume/sampler_kind.hpp"

namespace ttc
{

// Draws n points from {x : A x <= b}: every point is the endpoint of a fresh
// run of `walk` (walkLength steps) started at the certified Chebyshev centre,
// on volesti's mt19937 stream seeded with `seed`.  With `round` (ttc's
// default), the walk runs in the rounded body y (roundPolytope: maximum
// inscribed ellipsoid -> unit ball) and each point is mapped back,
// x = shift + T y; affine maps preserve uniformity.  Unrounded, every walk
// steps in units of the inner radius, so on a skinny body (R/r >> sqrt n) the
// samples never leave the centre (long box [0,2]x[0,1000]: 100% of samples in
// the last tenth), and a tiny body far from the origin loses its slacks to
// cancellation (samples outside P).  ballSeconds receives the time spent on
// the inner ball and the rounding.
std::vector<Eigen::VectorXd> sampleRestartWalk(const Eigen::MatrixXd& A,
                                               const Eigen::VectorXd& b, long long n,
                                               SamplerWalk walk, unsigned walkLength,
                                               unsigned seed, double& ballSeconds,
                                               bool round = true);

// GenerateSamples: inflate facets, independently restart the selected walk,
// stochastically round to 10^-b Z^n, and reject against the original body.
// Walk length and Volesti convergence remain practical, uncertified oracles.
template <typename Coord, typename GenFn>
std::vector<Eigen::Matrix<Coord, Eigen::Dynamic, 1>> sampleLattice(
    const Eigen::MatrixXd& A, const Eigen::VectorXd& b, long long n,
    int precision, double epsilon, double delta, unsigned seed,
    GenFn&& generate, double& samplingTime)
{
  using Vec = Eigen::Matrix<Coord, Eigen::Dynamic, 1>;
  using std::floor;
  using std::pow;
  const Coord scale = pow(Coord(10), precision);
  const double spacing = std::pow(10.0, -precision);
  if (spacing == 0.0)
    throw std::runtime_error("lattice spacing is below double geometry resolution");
  const double kappa = std::sqrt(2.0 * std::log(16.0 / epsilon))
                       + std::sqrt(2.0 * std::log(double(A.rows())));
  Eigen::VectorXd inflated = b;
  for (Eigen::Index i = 0; i < A.rows(); ++i)
    inflated(i) += spacing * kappa * A.row(i).norm();
  const auto Am = A.template cast<Coord>().eval();
  const Vec bm = b.template cast<Coord>();
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::vector<Vec> out;
  out.reserve(static_cast<std::size_t>(n));
  // Under the paper's oracle assumptions, rejection probability <= epsilon/3.
  // Budget a failed sampling batch within delta; report failure, never fewer points.
  const int maxBatches = std::max(1, static_cast<int>(std::ceil(
      (std::log(delta) - std::log(double(std::max(1LL, n)))) / std::log(epsilon / 3.0))));
  const auto begun = std::chrono::steady_clock::now();
  for (int batch = 0; out.size() < static_cast<std::size_t>(n); ++batch)
  {
    if (batch >= maxBatches)
      throw std::runtime_error("lattice rejection sampler exceeded its retry budget");
    double ignoredTime = 0.0;
    auto points = generate(A, inflated, n - static_cast<long long>(out.size()), ignoredTime);
    for (const auto& point : points)
    {
      Vec k(point.size());  // integer lattice coordinates: the point is k / 10^b
      for (Eigen::Index j = 0; j < point.size(); ++j)
      {
        const Coord raw = Coord(point(j)) * scale;
        Coord base = floor(raw);
        const Coord fraction = raw - base;
        if (unit(rng) < static_cast<double>(fraction)) base += 1;
        k(j) = base;
      }
      // Membership in lattice units, A k <= 10^b b.  Testing A (k / 10^b) <= b
      // decided lattice points that lie exactly on a facet by the rounding
      // error of k / 10^b (0.010 + 2*0.171 = 0.35200000000000004 > 0.352), so
      // such points were never produced.  In GMP both sides are exact; in
      // double allow a few ulps of the terms' magnitude.
      const Vec lhs = Am * k;
      const Vec rhs = bm * scale;
      bool inside;
      if constexpr (std::is_same<Coord, double>::value)
      {
        const Vec mag = Am.cwiseAbs() * k.cwiseAbs() + rhs.cwiseAbs();
        constexpr double tol = 8.0 * std::numeric_limits<double>::epsilon();
        inside = ((lhs - rhs).array() <= tol * mag.array()).all();
      }
      else
      {
        inside = ((lhs - rhs).array() <= Coord(0)).all();
      }
      if (inside) out.push_back(k / scale);
    }
  }
  samplingTime = std::chrono::duration<double>(std::chrono::steady_clock::now() - begun).count();
  return out;
}

}  // namespace ttc
