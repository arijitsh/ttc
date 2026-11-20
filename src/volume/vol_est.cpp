#include "volume/vol_est.hpp"

#include <algorithm>
#include <stdexcept>

#include "logger.hpp"
#include "volume/polytope_volume.hpp"
#include "volume/seed_streams.hpp"

namespace ttc
{

VolumeBounds ballVolumeBounds(double innerRadius, double outerRadius, std::size_t n,
                              double epsPrime)
{
  const double dim = double(n);
  return {(1.0 - epsPrime) * std::pow(2.0 * innerRadius / std::sqrt(dim), dim),
          (1.0 + epsPrime) * std::pow(2.0 * outerRadius, dim)};
}

VolEst::VolEst(VolumeOracle oracle, unsigned walkLength, bool round, double epsPrime, unsigned seed)
    : d_oracle(oracle), d_walkLength(walkLength), d_round(round), d_epsPrime(epsPrime), d_seed(seed)
{
}

double VolEst::operator()(const PreparedPolytope& P, std::size_t index)
{
  const std::size_t n = static_cast<std::size_t>(P.A.cols());
  const double start = Log.elapsed();
  // Round first (Cousins-Vempala, Fig. 2 step 1); vol(P) = det vol(rounded P)
  // is folded into RoundedPolytope.  P keeps its own radii for the clipping.
  const RoundedPolytope body(P.A, P.b, P.ball, d_round);
  if (body.rounded())
    Log(2) << "rounded polytope: |det T| = " << body.det() << ", inner radius " << P.ball.radius
           << " -> " << body.innerRadius() << std::endl;
  else if (d_round)
    Log(1) << "WARNING: polytope rounding skipped (" << body.roundingFailure()
           << "); the volume estimate may be poor if the polytope is skinny" << std::endl;

  constexpr unsigned kAttempts = 3;
  double volume = runUntilPositive(
      [&](unsigned attempt) {
        return body.volumeWithSeed(d_oracle, d_epsPrime, d_walkLength,
                                   oracleSeed(d_seed, index + 1, attempt));
      },
      kAttempts,
      [&](double failed) {
        Log(1) << "WARNING: volume oracle returned " << failed << " (inner radius "
               << P.ball.radius << ", outer radius " << P.outerRadius
               << "); retrying on a fresh random stream" << std::endl;
      });
  const double elapsed = Log.elapsed() - start;
  d_seconds += elapsed;
  Log(2) << "polytope volume " << volume << " (" << P.A.rows() << " facets, was " << P.rawFacets
         << ", dim " << n << ") in " << elapsed << "s" << std::endl;
  if (!std::isfinite(volume) || volume <= 0.0)
    throw std::runtime_error("volume oracle failed on a positive-dimensional polytope");

  const VolumeBounds bounds = ballVolumeBounds(P.ball.radius, P.outerRadius, n, d_epsPrime);
  if (std::isnan(bounds.lower) || std::isnan(bounds.upper) || bounds.lower > bounds.upper)
    throw std::runtime_error("invalid polytope volume clipping bounds");
  volume = std::clamp(volume, bounds.lower, bounds.upper);
  if (!std::isfinite(volume) || volume <= 0.0)
    throw std::runtime_error("clipped polytope volume exceeds double range");
  return volume;
}

}  // namespace ttc
