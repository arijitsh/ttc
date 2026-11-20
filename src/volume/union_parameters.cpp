#include "volume/union_parameters.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ttc
{

double lraSampleThreshold(double epsilon, double delta, std::size_t m)
{
  const double epsPrime = epsilon / 12.0;
  const double mm = std::max(1.0, static_cast<double>(m));
  return std::max(24.0 * (std::log(24.0) - std::log(delta)) / ((1.0 - epsPrime) * epsPrime * epsPrime),
                  15.0 * (std::log(6.0) - std::log(delta) + std::log(mm)));
}

UnionParameters unionParameters(double epsilon, double delta, std::size_t m, std::size_t n, int b)
{
  UnionParameters params;
  params.epsPrime = epsilon / 12.0;
  params.deltaPrime = delta / (4.0 * std::max(1.0, static_cast<double>(m)));
  params.thresh = lraSampleThreshold(epsilon, delta, m);
  if (!std::isfinite(params.thresh)
      || params.thresh >= double(std::numeric_limits<long long>::max() / 2)
      || params.deltaPrime <= 0.0)
    throw std::invalid_argument("epsilon/delta require an unsupported sample threshold");
  params.log2c = std::log2(params.thresh) + double(b) * double(n) * std::log2(10.0);
  return params;
}

}  // namespace ttc
