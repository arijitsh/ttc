#pragma once

// The single-polytope volume oracle selectable with --vol-oracle.  Kept free
// of volesti includes so lra_volume.hpp / main.cpp stay light; the oracle
// itself lives in polytope_volume_oracle.hpp.

#include <cstddef>

namespace ttc
{

// Volume oracle run on each (rounded) polytope.  Selected with --vol-oracle.
//   CoolingGaussiansBall (default, 'cg-ball'): Gaussian cooling with the
//     Gaussian ball walk -- the Cousins-Vempala annealing schedule.
//   CoolingGaussiansCdhr ('cg-cdhr'): same schedule, coordinate hit-and-run
//     (the walk of Cousins-Vempala's reference implementation).
//   CoolingBallsBilliard ('cb-billiard'): volesti's cooling balls with the
//     billiard walk; much faster, somewhat noisier.
enum class VolumeOracle
{
  CoolingGaussiansBall,
  CoolingGaussiansCdhr,
  CoolingBallsBilliard
};

inline const char* volumeOracleName(VolumeOracle oracle)
{
  switch (oracle)
  {
    case VolumeOracle::CoolingGaussiansCdhr: return "cg-cdhr";
    case VolumeOracle::CoolingBallsBilliard: return "cb-billiard";
    case VolumeOracle::CoolingGaussiansBall:
    default: return "cg-ball";
  }
}

// Default walk length per oracle (the settings validated in
// tests/volesti/test_alt_oracles.cpp): the ball walk needs ~n^2 steps per
// sample, coordinate hit-and-run volesti's default 10 + n/10, and one billiard
// trajectory per sample suffices.
inline unsigned defaultOracleWalkLength(VolumeOracle oracle, std::size_t dimension)
{
  switch (oracle)
  {
    case VolumeOracle::CoolingGaussiansCdhr:
      return static_cast<unsigned>(10 + dimension / 10);
    case VolumeOracle::CoolingBallsBilliard: return 1u;
    case VolumeOracle::CoolingGaussiansBall:
    default: return static_cast<unsigned>(10 + dimension * dimension);
  }
}

}  // namespace ttc
