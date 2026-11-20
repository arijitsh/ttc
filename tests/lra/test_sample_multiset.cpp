// The multiset X and the probability p of Algorithm 1
// (src/volume/sample_multiset.hpp, union_estimator.hpp):
//   lines 9-11  remove the points of X inside P_i
//   lines 13/17 keep each point of X with probability 1/2
//   line 20     append the new samples
//   line 22     |X| / (p c), with p tracked as log2(p c)

#include <cmath>
#include <random>
#include <vector>

#include "lra_test_util.hpp"
#include "volume/mp_float.hpp"
#include "volume/sample_multiset.hpp"
#include "volume/union_estimator.hpp"

namespace
{

template <typename Coord>
using Vec = typename ttc::SampleMultiset<Coord>::Vec;

template <typename Coord>
Vec<Coord> point(double x, double y)
{
  Vec<Coord> v(2);
  v << Coord(x), Coord(y);
  return v;
}

template <typename Coord>
void removeInsideKeepsExactlyTheOutsidePoints()
{
  ttc::SampleMultiset<Coord> X;
  // Inside, on the boundary (counts as inside), outside, and a duplicate.
  X.append(std::vector<Vec<Coord>>{point<Coord>(0.5, 0.5), point<Coord>(1.0, 0.25), point<Coord>(0.0, 0.0),
                                   point<Coord>(1.5, 0.5), point<Coord>(0.5, -0.1), point<Coord>(0.5, 0.5)});
  REQUIRE(X.size() == 6);
  const auto B = lra_test::box({0, 0}, {1, 1});
  CHECK(X.removeInside(B.A, B.b) == 4);
  REQUIRE(X.size() == 2);
  for (const auto& s : X.points()) CHECK(((B.A.template cast<Coord>() * s - B.b.template cast<Coord>()).array() > Coord(0)).any());
  CHECK(X.removeInside(Eigen::MatrixXd(0, 2), Eigen::VectorXd(0)) == 0);  // no rows: nothing removed
}

}  // namespace

TEST_CASE("X, lines 9-11: removeInside drops exactly the points of P_i (boundary included)")
{
  SUBCASE("double") { removeInsideKeepsExactlyTheOutsidePoints<double>(); }
  SUBCASE("gmp")
  {
    ttc::MpFloat::default_precision(50);
    removeInsideKeepsExactlyTheOutsidePoints<ttc::MpFloat>();
  }
}

TEST_CASE("X, lines 13 and 17: halve keeps each point with probability 1/2, one draw per point")
{
  ttc::SampleMultiset<double> X;
  const int n = 20000;
  X.append(std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Zero(1)));
  std::mt19937 rng(11), replay(11);
  const std::size_t removed = X.halve(rng);
  CHECK(removed + X.size() == std::size_t(n));
  CHECK(std::abs(double(X.size()) - n / 2.0) < 4.0 * std::sqrt(n / 4.0));  // 4 sigma
  // Exactly one uniform per point: the streams line up again afterwards.
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (int i = 0; i < n; ++i) unit(replay);
  CHECK(rng() == replay());

  ttc::SampleMultiset<double> empty;
  CHECK(empty.halve(rng) == 0);
}

TEST_CASE("X, line 22: |X| / (p c) from log2(p c)")
{
  // c = 1000, three halvings: p = 1/8, so the estimate is |X| 8 / 1000.
  ttc::ScaledProbability p(std::log2(1000.0));
  for (int i = 0; i < 3; ++i) p.halve();
  CHECK(p.estimate(10) == doctest::Approx(0.08));
  CHECK(p.estimate(0) == 0.0);
  // log2(t p) for t = c vol: vol = 2 -> t p = 2000 / 8.
  CHECK(std::exp2(p.log2tp(2.0)) == doctest::Approx(250.0));

  // Huge c (b n = 300) and tiny volumes stay finite.
  ttc::ScaledProbability big(std::log2(27699.0) + 300 * std::log2(10.0));
  CHECK(std::isfinite(big.log2tp(1e-290)));
  CHECK(std::abs(big.estimate(27699) / 1e-300 - 1.0) < 1e-9);  // relative: Approx is absolute near 0
}
