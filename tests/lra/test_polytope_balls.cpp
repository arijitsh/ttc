// The balls B(z, gamma) ⊆ P ⊆ B(z, R) (src/volume/polytope_balls.hpp):
// gamma is MaxInscribedBall of GetPrecision (Alg. 2, line 7), and both radii
// bound the VolEst clipping of line 8.

#include <cmath>
#include <stdexcept>

#include "lra_test_util.hpp"
#include "volume/polytope_balls.hpp"

TEST_CASE("balls: the Chebyshev ball of a box and of a simplex")
{
  const auto B = lra_test::box({0, 0, 0}, {2, 6, 4});
  const auto ball = ttc::certifiedInnerBall(B.A, B.b);
  CHECK(ball.radius == doctest::Approx(1.0).epsilon(1e-6));
  CHECK(((B.A * ball.center - B.b).array() <= -ball.radius + 1e-6).all());

  // The standard n-simplex has inradius 1 / (n + sqrt n).
  for (int n : {2, 5, 10})
  {
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    lra_test::simplex(n, A, b);
    CAPTURE(n);
    CHECK(ttc::certifiedInnerBall(A, b).radius == doctest::Approx(1.0 / (n + std::sqrt(double(n)))).epsilon(1e-6));
  }
}

TEST_CASE("balls: empty and flat polytopes have no inner ball")
{
  Eigen::MatrixXd A(4, 2);
  Eigen::VectorXd b(4);
  A << 1, 0, -1, 0, 0, 1, 0, -1;
  b << 1, -1, 1, 0;  // x = 1, 0 <= y <= 1
  CHECK(ttc::certifiedInnerBall(A, b).radius <= 1e-9);
  b << 0, -1, 1, 0;  // x <= 0 and x >= 1
  CHECK(ttc::certifiedInnerBall(A, b).radius <= 0.0);
}

TEST_CASE("balls: the outer radius contains every vertex and is not loose")
{
  const auto B = lra_test::box({-1, 3}, {5, 4});
  const auto ball = ttc::certifiedInnerBall(B.A, B.b);
  const double R = ttc::outerRadius(B.A, B.b, ball.center);
  double farthest = 0.0;
  for (double x : {-1.0, 5.0})
    for (double y : {3.0, 4.0}) farthest = std::max(farthest, (Eigen::Vector2d(x, y) - ball.center).norm());
  CHECK(R >= farthest - 1e-9);
  CHECK(R <= std::sqrt(2.0) * farthest + 1e-9);  // a box bound is within sqrt(n) of the true one
}

TEST_CASE("balls: an unbounded polytope has no outer radius")
{
  Eigen::MatrixXd A(2, 2);
  Eigen::VectorXd b(2);
  A << -1, 0, 0, -1;  // the positive quadrant
  b << 0, 0;
  CHECK_THROWS_AS(ttc::outerRadius(A, b, Eigen::Vector2d(1, 1)), std::runtime_error);
}
