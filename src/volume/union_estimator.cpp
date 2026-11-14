#include "volume/union_estimator.hpp"

#include <limits>
#include <stdexcept>

namespace ttc
{

// GetPrecision in the revised manuscript: eta = epsilon/96. Work in logs
// so gamma/minRadius need not be representable as a double.
int precisionFromCubes(std::size_t dimension, std::size_t facetCount,
                       std::size_t numCubes, double minRadius, double epsilon)
{
  const double eta = epsilon / 96.0;
  const double logGamma = std::log(64.0) + std::log(double(numCubes))
      + std::log(double(dimension)) - std::log(eta)
      + 0.5 * std::log(std::log(4.0) + std::log(double(facetCount)) - std::log(eta));
  const double prec = std::max(0.0, std::ceil((logGamma - std::log(minRadius))
                                           / std::log(10.0)));
  if (!std::isfinite(prec) || prec > std::numeric_limits<int>::max() - 32)
  {
    throw std::runtime_error("lattice precision exceeds supported range");
  }
  return static_cast<int>(prec);
}

double lraSampleThreshold(double epsilon, double delta, std::size_t numPolytopes)
{
  const double mvcEps = epsilon / 12.0;
  const double m = std::max(1.0, static_cast<double>(numPolytopes));
  return std::max(
      24.0 * (std::log(24.0) - std::log(delta)) / ((1.0 - mvcEps) * mvcEps * mvcEps),
      15.0 * (std::log(6.0) - std::log(delta) + std::log(m)));
}

}  // namespace ttc
