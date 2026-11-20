// Line 3 of Algorithm 1, b <- GetPrecision(F_bd, M, eps/96)  (Algorithm 2,
// src/volume/get_precision.hpp).

#include <cmath>
#include <vector>

#include "lra_test_util.hpp"
#include "volume/get_precision.hpp"

namespace
{

// Algorithm 2 line by line, in long double.
int algorithm2(std::size_t n, std::size_t m, const std::vector<long>& facets,
               const std::vector<double>& gammaHat, double eta)
{
  long double maxratio = 1.0L;
  long double r = 0.0L;
  for (long f : facets) r += f;
  for (double g : gammaHat)
  {
    const long double gamma = 64.0L * m * n / eta * std::sqrt(std::log(8.0L * r / eta));
    maxratio = std::max(maxratio, gamma / g);
  }
  return static_cast<int>(std::ceil(std::log10(maxratio)));
}

ttc::PreparedPolytope polytopeWith(long facets, double innerRadius)
{
  ttc::PreparedPolytope P;
  P.A = Eigen::MatrixXd::Zero(facets, 2);
  P.b = Eigen::VectorXd::Zero(facets);
  P.ball.center = Eigen::Vector2d::Zero();
  P.ball.radius = innerRadius;
  P.outerRadius = 1.0;
  return P;
}

}  // namespace

TEST_CASE("GetPrecision: equals Algorithm 2 on a grid of inputs")
{
  for (std::size_t n : {1u, 3u, 10u, 40u})
    for (std::size_t m : {1u, 4u, 100u})
      for (double eps : {0.1, 0.8, 1.0})
        for (double gamma : {1e-6, 0.37, 1.0, 1e3})
        {
          const double eta = eps / 96.0;
          const std::vector<long> facets(m, 2 * long(n));
          const std::vector<double> radii(m, gamma);
          CAPTURE(n);
          CAPTURE(m);
          CAPTURE(eps);
          CAPTURE(gamma);
          CHECK(ttc::latticePrecision(n, 2 * n * m, m, gamma, eta) == algorithm2(n, m, facets, radii, eta));
        }
}

TEST_CASE("GetPrecision: r is the total facet count and gamma_hat the smallest inner radius")
{
  const std::vector<ttc::PreparedPolytope> P = {polytopeWith(4, 2.0), polytopeWith(7, 0.05),
                                                polytopeWith(3, 9.0)};
  const double eta = 0.8 / 96.0;
  const std::size_t m = 5;  // m counts every cube, kept or not
  CHECK(ttc::getPrecision(P, 2, m, eta) == ttc::latticePrecision(2, 14, m, 0.05, eta));
  CHECK(ttc::getPrecision(P, 2, m, eta) == algorithm2(2, m, {4, 7, 3}, {2.0, 0.05, 9.0}, eta));
  CHECK(ttc::getPrecision({}, 2, m, eta) == 0);
}

TEST_CASE("GetPrecision: b = 0 for a large ball, and +1 per factor 10 smaller")
{
  const double eta = 0.8 / 96.0;
  CHECK(ttc::latticePrecision(2, 4, 1, 1e12, eta) == 0);
  const int b = ttc::latticePrecision(3, 6, 2, 0.3, eta);
  REQUIRE(b > 0);
  CHECK(ttc::latticePrecision(3, 6, 2, 0.03, eta) == b + 1);
  CHECK(ttc::latticePrecision(3, 6, 2, 3e-7, eta) == b + 6);
  // Smaller eta (higher accuracy) never lowers b.
  CHECK(ttc::latticePrecision(3, 6, 2, 0.3, eta / 10) >= b);
}
