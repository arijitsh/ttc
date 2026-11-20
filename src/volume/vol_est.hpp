#pragma once

// Line 8 of Algorithm 1:  t <- Thresh 10^{bn} VolEst(P_i, eps', delta').
// This is the VolEst factor; the loop (union_estimator.hpp) applies c.
//
// VolEst runs a practical Cousins-Vempala oracle (polytope_volume.hpp) on the
// rounded polytope.  It takes no delta (VolEsti has no confidence parameter),
// so delta' is not used.  A run that ends non-finite or non-positive is a
// detected failure and is retried on a fresh stream.  The estimate is then
// clipped to the ball bounds
//     [(1 - eps') (2 gamma_i / sqrt n)^n,  (1 + eps') (2 R_i)^n],
// which hold for every polytope with B(gamma_i) ⊆ P_i ⊆ B(R_i) and keep the
// number of halvings in lines 12-14 polynomial.

#include <cmath>
#include <cstddef>

#include "volume/prepared_polytope.hpp"
#include "volume/volume_oracle_kind.hpp"

namespace ttc
{

struct VolumeBounds
{
  double lower = 0.0;
  double upper = 0.0;
};

// The clipping interval for inner radius gamma and outer radius R in R^n.
VolumeBounds ballVolumeBounds(double innerRadius, double outerRadius, std::size_t n,
                              double epsPrime);

// Runs oracle(0), oracle(1), ... until one returns a finite positive volume,
// at most `attempts` times, calling onRetry(failedValue) before each rerun.
// Returns the last value.
template <typename Oracle, typename OnRetry>
double runUntilPositive(Oracle&& oracle, unsigned attempts, OnRetry&& onRetry)
{
  double volume = oracle(0u);
  for (unsigned attempt = 1; attempt < attempts && !(std::isfinite(volume) && volume > 0.0);
       ++attempt)
  {
    onRetry(volume);
    volume = oracle(attempt);
  }
  return volume;
}

class VolEst
{
 public:
  VolEst(VolumeOracle oracle, unsigned walkLength, bool round, double epsPrime, unsigned seed);

  // The clipped estimate of vol(P); `index` is P's 0-based position among the
  // kept polytopes and selects its random streams.  Throws when every oracle
  // run fails or the bounds are not representable.
  double operator()(const PreparedPolytope& P, std::size_t index);

  double seconds() const { return d_seconds; }  // total time in the oracle

 private:
  VolumeOracle d_oracle;
  unsigned d_walkLength;
  bool d_round;
  double d_epsPrime;
  unsigned d_seed;
  double d_seconds = 0.0;
};

}  // namespace ttc
