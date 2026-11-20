#pragma once

// Line 19 of Algorithm 1:  S <- GenerateSamples(P, N, b, eps', delta').
//
// N lattice points of P ∩ 10^-b Z^n, each within a (1 ± eps') factor of
// uniform (proof appendix, cor:linf, after Kannan-Vempala 1997):
//
//   1  P' <- P with every facet j pushed out by kappa 10^-b ||A_j||,
//          kappa = sqrt(2 ln(16/eps')) + sqrt(2 ln r)          (r = facets of P)
//   2  repeat
//   3    p <- a point of P' from the continuous sampler
//   4    x <- p rounded at random to 10^-b Z^n, coordinatewise with E[x] = p
//   5  until x ∈ P
//
// Every sample is its own run of 2-5.  The continuous sampler `proposals`
// (polytope_sampler.hpp: a walk restarted at the centre of P') is asked for
// all missing points at once, in batches.  Under the oracle assumptions a
// proposal is rejected with probability <= eps'/3, so the batch budget makes
// running out a probability-<= delta' event; the sampler then throws rather
// than return fewer than N points.  The walk length is a practical choice, not
// a certified mixing time.

#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>

namespace ttc
{
namespace lattice
{

// kappa of step 1, for accuracy eps' and r facets.
inline double facetInflation(double epsPrime, long facets)
{
  return std::sqrt(2.0 * std::log(16.0 / epsPrime)) + std::sqrt(2.0 * std::log(double(facets)));
}

// The right-hand side of P' (step 1): b_j + distance ||A_j||.
inline Eigen::VectorXd inflateFacets(const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
                                     double distance)
{
  Eigen::VectorXd inflated = b;
  for (Eigen::Index j = 0; j < A.rows(); ++j) inflated(j) += distance * A.row(j).norm();
  return inflated;
}

// Batches allowed before giving up: the smallest B with (eps'/3)^B <= delta'/N.
inline int batchBudget(double epsPrime, double deltaPrime, long long N)
{
  return std::max(1, static_cast<int>(std::ceil((std::log(deltaPrime) - std::log(double(std::max(1LL, N))))
                                                / std::log(epsPrime / 3.0))));
}

// Step 4 in lattice units: k with k_j = floor(p_j 10^b) + Bernoulli(frac(p_j 10^b)),
// so that k / 10^b is the rounded point.  One uniform draw per coordinate.
template <typename Coord, typename Point, typename URNG>
Eigen::Matrix<Coord, Eigen::Dynamic, 1> randomizedRound(const Point& p, const Coord& scale,
                                                        URNG& rng)
{
  using std::floor;
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  Eigen::Matrix<Coord, Eigen::Dynamic, 1> k(p.size());
  for (Eigen::Index j = 0; j < p.size(); ++j)
  {
    const Coord raw = Coord(p(j)) * scale;
    Coord base = floor(raw);
    const Coord fraction = raw - base;
    if (unit(rng) < static_cast<double>(fraction)) base += 1;
    k(j) = base;
  }
  return k;
}

// Step 5 in lattice units, A k <= 10^b b (scaledB = 10^b b).  Testing
// A (k / 10^b) <= b instead decided lattice points that lie exactly on a
// facet by the rounding error of k / 10^b (0.010 + 2*0.171 =
// 0.35200000000000004 > 0.352), so such points were never produced.  In GMP
// both sides are exact; in double a few ulps of the terms' magnitude are allowed.
template <typename Coord>
bool latticePointInside(const Eigen::Matrix<Coord, Eigen::Dynamic, Eigen::Dynamic>& A,
                        const Eigen::Matrix<Coord, Eigen::Dynamic, 1>& scaledB,
                        const Eigen::Matrix<Coord, Eigen::Dynamic, 1>& k)
{
  using Vec = Eigen::Matrix<Coord, Eigen::Dynamic, 1>;
  const Vec lhs = A * k;
  if constexpr (std::is_same<Coord, double>::value)
  {
    const Vec magnitude = A.cwiseAbs() * k.cwiseAbs() + scaledB.cwiseAbs();
    constexpr double tol = 8.0 * std::numeric_limits<double>::epsilon();
    return ((lhs - scaledB).array() <= tol * magnitude.array()).all();
  }
  else
  {
    return ((lhs - scaledB).array() <= Coord(0)).all();
  }
}

}  // namespace lattice

// GenerateSamples(P, N, b, eps', delta') with P = {x : A x <= rhs}.
//   Coord      arithmetic of the rounding, the membership test and the result
//   proposals  (A, b', count, ignoredSeconds) -> `count` points of {A x <= b'}
//   seed       the stream of the randomized rounding
// samplingSeconds receives the time of the whole call.
template <typename Coord, typename Proposals>
std::vector<Eigen::Matrix<Coord, Eigen::Dynamic, 1>> generateSamples(
    const Eigen::MatrixXd& A, const Eigen::VectorXd& rhs, long long N, int b, double epsPrime,
    double deltaPrime, unsigned seed, Proposals&& proposals, double& samplingSeconds)
{
  using Vec = Eigen::Matrix<Coord, Eigen::Dynamic, 1>;
  using std::pow;
  const Coord scale = pow(Coord(10), b);  // 10^b
  const double spacing = std::pow(10.0, -b);
  if (spacing == 0.0) throw std::runtime_error("lattice spacing is below double geometry resolution");

  const double kappa = lattice::facetInflation(epsPrime, A.rows());
  const Eigen::VectorXd inflatedB = lattice::inflateFacets(A, rhs, spacing * kappa);  // step 1
  const auto Am = A.template cast<Coord>().eval();
  const Vec scaledB = Vec(rhs.template cast<Coord>()) * scale;
  std::mt19937 rng(seed);
  const int maxBatches = lattice::batchBudget(epsPrime, deltaPrime, N);

  std::vector<Vec> S;
  S.reserve(static_cast<std::size_t>(N));
  const auto begun = std::chrono::steady_clock::now();
  for (int batch = 0; S.size() < static_cast<std::size_t>(N); ++batch)  // steps 2-5
  {
    if (batch >= maxBatches)
      throw std::runtime_error("lattice rejection sampler exceeded its retry budget");
    double ignoredSeconds = 0.0;
    const auto points = proposals(A, inflatedB, N - static_cast<long long>(S.size()), ignoredSeconds);
    for (const auto& p : points)  // step 3
    {
      const Vec k = lattice::randomizedRound<Coord>(p, scale, rng);  // step 4
      if (lattice::latticePointInside<Coord>(Am, scaledB, k)) S.push_back(k / scale);  // step 5
    }
  }
  samplingSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begun).count();
  return S;
}

}  // namespace ttc
