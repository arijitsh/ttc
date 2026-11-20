// Line 19 of Algorithm 1, S <- GenerateSamples(P_i, N_i, b, eps', delta')
// (src/volume/generate_samples.hpp): inflate P to P', propose from P',
// round at random to 10^-b Z^n, reject outside P.

#include <cmath>
#include <map>
#include <random>
#include <stdexcept>
#include <vector>

#include "lra_test_util.hpp"
#include "volume/generate_samples.hpp"
#include "volume/mp_float.hpp"
#include "volume/polytope_sampler.hpp"

namespace lat = ttc::lattice;

namespace
{

const double kEpsPrime = 0.8 / 12.0;
const double kDeltaPrime = 0.05;

// True if every coordinate of x is an integer multiple of 10^-b.
bool onLattice(const Eigen::VectorXd& x, int b)
{
  for (Eigen::Index j = 0; j < x.size(); ++j)
  {
    const double k = x(j) * std::pow(10.0, b);
    if (std::abs(k - std::round(k)) > 1e-6) return false;
  }
  return true;
}

}  // namespace

TEST_CASE("GenerateSamples, step 1: kappa and the inflated facets of P'")
{
  CHECK(lat::facetInflation(kEpsPrime, 6)
        == doctest::Approx(std::sqrt(2 * std::log(16 / kEpsPrime)) + std::sqrt(2 * std::log(6.0))));
  Eigen::MatrixXd A(2, 2);
  A << 3, 4, -1, 0;
  Eigen::VectorXd b(2);
  b << 1, 2;
  const Eigen::VectorXd inflated = lat::inflateFacets(A, b, 0.1);
  CHECK(inflated(0) == doctest::Approx(1 + 0.1 * 5));  // ||(3,4)|| = 5
  CHECK(inflated(1) == doctest::Approx(2 + 0.1 * 1));
}

TEST_CASE("GenerateSamples, budget: the fewest batches B with (eps'/3)^B <= delta'/N")
{
  for (long long N : {1LL, 10LL, 1000LL, 100000LL})
  {
    const int B = lat::batchBudget(kEpsPrime, kDeltaPrime, N);
    CAPTURE(N);
    CHECK(std::pow(kEpsPrime / 3, B) <= kDeltaPrime / N * (1 + 1e-12));
    if (B > 1) CHECK(std::pow(kEpsPrime / 3, B - 1) > kDeltaPrime / N);
  }
}

TEST_CASE("GenerateSamples, step 4: randomized rounding goes to a neighbour, unbiased")
{
  std::mt19937 rng(5);
  const double scale = 100.0;  // b = 2
  Eigen::VectorXd p(2);
  p << 0.12345, -3.0;  // -3.0 is a lattice point and must stay put
  double sum = 0.0;
  const int trials = 40000;
  for (int t = 0; t < trials; ++t)
  {
    const Eigen::VectorXd k = lat::randomizedRound<double>(p, scale, rng);
    CHECK((k(0) == 12.0 || k(0) == 13.0));
    CHECK(k(1) == -300.0);
    sum += k(0);
  }
  // E[k] = 12.345; the sd of the mean is sqrt(0.345 * 0.655 / trials).
  CHECK(std::abs(sum / trials - 12.345) < 4 * std::sqrt(0.345 * 0.655 / trials));
}

TEST_CASE("GenerateSamples, step 5: lattice points on a facet are inside, in double and GMP")
{
  // 0.010 + 2 * 0.171 = 0.352 exactly, but not in double arithmetic on x itself.
  Eigen::MatrixXd A(1, 2);
  A << 1, 2;
  const double scale = 1000.0;
  Eigen::VectorXd k(2);
  k << 10, 171;
  Eigen::VectorXd scaledB(1);
  scaledB << 0.352 * scale;
  CHECK(lat::latticePointInside<double>(A, scaledB, k));
  k(1) = 172;
  CHECK_FALSE(lat::latticePointInside<double>(A, scaledB, k));

  ttc::MpFloat::default_precision(50);
  using MpMat = Eigen::Matrix<ttc::MpFloat, Eigen::Dynamic, Eigen::Dynamic>;
  const MpMat Am = A.cast<ttc::MpFloat>();
  ttc::MpVector km(2), bm(1);
  km << 10, 171;
  bm << ttc::MpFloat(352);
  CHECK(lat::latticePointInside<ttc::MpFloat>(Am, bm, km));
  km(1) = 172;
  CHECK_FALSE(lat::latticePointInside<ttc::MpFloat>(Am, bm, km));
}

TEST_CASE("GenerateSamples: N lattice points of P, uniform over P ∩ 10^-b Z^n")
{
  // The analysis needs P' barely larger than P (kappa 10^-b << inner radius),
  // so use [0,1] at b = 3: 1001 lattice points and P' = [-0.005, 1.005].
  // With exact uniform proposals each point appears about N / 1001 times.
  const auto B = lra_test::box({0}, {1});
  const long long perPoint = 200;
  const long long N = 1001 * perPoint;
  lra_test::UniformBoxSampler proposals;
  double seconds = 0.0;
  const auto S = ttc::generateSamples<double>(B.A, B.b, N, 3, kEpsPrime, kDeltaPrime, 9u, proposals, seconds);
  REQUIRE(S.size() == std::size_t(N));
  std::map<long, int> count;
  for (const auto& x : S)
  {
    REQUIRE(onLattice(x, 3));
    REQUIRE(B.contains(x));
    ++count[std::lround(x(0) * 1000)];
  }
  CHECK(count.size() == 1001);
  const double sd = std::sqrt(double(perPoint));
  for (const auto& [k, c] : count)
  {
    CAPTURE(k);
    CHECK(std::abs(c - perPoint) < 5.5 * sd);
  }
  // The end points get the same mass as interior points (the rounding tent fits in P').
  CHECK(std::abs(count[0] + count[1000] - 2 * perPoint) < 5.5 * std::sqrt(2.0) * sd);

  // In 2D at b = 3 every sample is a lattice point of the box.
  const auto B2 = lra_test::box({0, 0}, {1, 1});
  for (const auto& x : ttc::generateSamples<double>(B2.A, B2.b, 5000, 3, kEpsPrime, kDeltaPrime, 4u, proposals, seconds))
  {
    CHECK(onLattice(x, 3));
    CHECK(B2.contains(x));
  }
}

TEST_CASE("GenerateSamples: N = 0, and a sampler that never lands in P exhausts its budget")
{
  const auto B = lra_test::box({0, 0}, {1, 1});
  int calls = 0;
  auto never = [&](const Eigen::MatrixXd&, const Eigen::VectorXd&, long long count, double&) {
    ++calls;
    return std::vector<Eigen::VectorXd>(std::size_t(count), Eigen::Vector2d(5, 5));
  };
  double seconds = 0.0;
  CHECK(ttc::generateSamples<double>(B.A, B.b, 0, 2, kEpsPrime, kDeltaPrime, 1u, never, seconds).empty());
  CHECK(calls == 0);
  CHECK_THROWS_AS(ttc::generateSamples<double>(B.A, B.b, 10, 2, kEpsPrime, kDeltaPrime, 1u, never, seconds),
                  std::runtime_error);
  CHECK(calls == lat::batchBudget(kEpsPrime, kDeltaPrime, 10));
}

TEST_CASE("GenerateSamples with ttc's restart walk: lattice points of a skinny triangle")
{
  // The real proposal sampler (billiard walk, rounded body), in double and GMP.
  Eigen::MatrixXd A(3, 2);
  A << -1, 0, 0, -1, 1, 50;  // x, y >= 0, x + 50 y <= 10
  Eigen::VectorXd b(3);
  b << 0, 0, 10;
  auto walk = [](const Eigen::MatrixXd& Ai, const Eigen::VectorXd& bi, long long count, double&) {
    double ballSeconds = 0.0;
    return ttc::sampleRestartWalk(Ai, bi, count, ttc::SamplerWalk::Billiard, 2, 77u, ballSeconds);
  };
  double seconds = 0.0;
  const auto S = ttc::generateSamples<double>(A, b, 500, 4, kEpsPrime, kDeltaPrime, 3u, walk, seconds);
  REQUIRE(S.size() == 500);
  for (const auto& x : S)
  {
    CHECK(onLattice(x, 4));
    CHECK(((A * x - b).array() <= 1e-12).all());
  }

  ttc::MpFloat::default_precision(40);
  const auto G = ttc::generateSamples<ttc::MpFloat>(A, b, 200, 4, kEpsPrime, kDeltaPrime, 3u, walk, seconds);
  REQUIRE(G.size() == 200);
  for (const auto& x : G) CHECK(((A.cast<ttc::MpFloat>() * x - b.cast<ttc::MpFloat>()).array() <= 0).all());
}
