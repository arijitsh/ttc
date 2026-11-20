// Line 8 of Algorithm 1, VolEst(P_i, eps', delta') with retries and clipping
// (src/volume/vol_est.hpp).

#include <cmath>
#include <limits>
#include <vector>

#include "lra_test_util.hpp"
#include "volume/polytope_balls.hpp"
#include "volume/vol_est.hpp"

namespace
{

ttc::PreparedPolytope prepare(const Eigen::MatrixXd& A, const Eigen::VectorXd& b)
{
  ttc::PreparedPolytope P;
  P.A = A;
  P.b = b;
  P.ball = ttc::certifiedInnerBall(A, b);
  P.outerRadius = ttc::outerRadius(A, b, P.ball.center);
  P.rawFacets = A.rows();
  return P;
}

ttc::VolEst defaultVolEst(std::size_t n, double eps = 0.8)
{
  return ttc::VolEst(ttc::kDefaultVolumeOracle, ttc::defaultOracleWalkLength(ttc::kDefaultVolumeOracle, n),
                     true, eps / 12.0, 123);
}

}  // namespace

TEST_CASE("VolEst: the clip bounds hold for boxes and simplices")
{
  const double epsPrime = 0.8 / 12.0;
  std::vector<std::pair<ttc::PreparedPolytope, double>> cases;  // (P, exact volume)
  for (const auto& B : {lra_test::box({0, 0}, {1, 1}), lra_test::box({-3, 0, 1}, {5, 0.1, 2}),
                        lra_test::box({0, 0, 0, 0, 0}, {1, 2, 3, 4, 5})})
    cases.emplace_back(prepare(B.A, B.b), B.volume());
  for (int n : {2, 4, 8})
  {
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    lra_test::simplex(n, A, b);
    cases.emplace_back(prepare(A, b), 1.0 / std::tgamma(n + 1.0));
  }
  for (const auto& [P, exact] : cases)
  {
    const auto bounds = ttc::ballVolumeBounds(P.ball.radius, P.outerRadius, P.A.cols(), epsPrime);
    CAPTURE(exact);
    CHECK(bounds.lower <= exact);
    CHECK(exact <= bounds.upper);
  }
}

TEST_CASE("VolEst: retries a failed oracle run on fresh attempts, at most `attempts` times")
{
  std::vector<unsigned> calls;
  std::vector<double> failures;
  auto onRetry = [&](double v) { failures.push_back(v); };

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::vector<double> answers = {nan, -1.0, 7.5};
  auto oracle = [&](unsigned attempt) {
    calls.push_back(attempt);
    return answers[attempt];
  };
  CHECK(ttc::runUntilPositive(oracle, 3, onRetry) == 7.5);
  CHECK(calls == std::vector<unsigned>{0, 1, 2});
  CHECK(failures.size() == 2);

  calls.clear();
  CHECK(ttc::runUntilPositive(oracle, 2, onRetry) == -1.0);  // gives up: returns the last value
  CHECK(calls == std::vector<unsigned>{0, 1});

  calls.clear();
  CHECK(ttc::runUntilPositive([&](unsigned a) { calls.push_back(a); return 2.0; }, 3, onRetry) == 2.0);
  CHECK(calls.size() == 1);
}

TEST_CASE("VolEst: estimates simple bodies within a few eps'")
{
  const auto cube = lra_test::box({0, 0, 0}, {1, 1, 1});
  const auto skinny = lra_test::box({0, 0}, {2, 1000});  // needs the rounding step
  Eigen::MatrixXd S;
  Eigen::VectorXd s;
  lra_test::simplex(4, S, s);
  struct Case { ttc::PreparedPolytope P; double exact; };
  const std::vector<Case> cases = {{prepare(cube.A, cube.b), 1.0},
                                   {prepare(skinny.A, skinny.b), 2000.0},
                                   {prepare(S, s), 1.0 / 24.0}};
  for (std::size_t i = 0; i < cases.size(); ++i)
  {
    auto volEst = defaultVolEst(cases[i].P.A.cols());
    const double v = volEst(cases[i].P, i);
    CAPTURE(i);
    CHECK(v == doctest::Approx(cases[i].exact).epsilon(3 * 0.8 / 12.0));
    CHECK(volEst.seconds() > 0.0);
  }
}

TEST_CASE("VolEst: same seed and index give the same estimate; another index another stream")
{
  const auto B = lra_test::box({0, 0, 0}, {1, 2, 3});
  const auto P = prepare(B.A, B.b);
  auto first = defaultVolEst(3);
  auto second = defaultVolEst(3);
  const double v0 = first(P, 0);
  CHECK(second(P, 0) == v0);
  CHECK(first(P, 1) != v0);
}
