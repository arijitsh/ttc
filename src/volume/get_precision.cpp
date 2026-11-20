#include "volume/get_precision.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ttc
{

int getPrecision(const std::vector<PreparedPolytope>& P, std::size_t n, std::size_t m, double eta)
{
  if (P.empty()) return 0;
  std::size_t r = 0;  // line 2
  double minInnerRadius = std::numeric_limits<double>::infinity();  // line 7, over i
  for (const PreparedPolytope& Pi : P)
  {
    r += static_cast<std::size_t>(Pi.A.rows());
    minInnerRadius = std::min(minInnerRadius, Pi.ball.radius);
  }
  return latticePrecision(n, r, m, minInnerRadius, eta);
}

int latticePrecision(std::size_t n, std::size_t r, std::size_t m, double minInnerRadius,
                     double eta)
{
  // ln gamma = ln(64 m n / eta) + ln(ln(8 r / eta)) / 2
  const double logGamma = std::log(64.0) + std::log(double(m)) + std::log(double(n)) - std::log(eta)
                          + 0.5 * std::log(std::log(8.0) + std::log(double(r)) - std::log(eta));
  const double b =
      std::max(0.0, std::ceil((logGamma - std::log(minInnerRadius)) / std::log(10.0)));
  if (!std::isfinite(b) || b > std::numeric_limits<int>::max() - 32)
    throw std::runtime_error("lattice precision exceeds supported range");
  return static_cast<int>(b);
}

}  // namespace ttc
