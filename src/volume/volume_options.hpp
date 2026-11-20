#pragma once

// Options of the LRA volume engine (unionvolume.hpp) and the settings it
// resolves before the main loop.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>

#include "volume/sampler_kind.hpp"
#include "volume/volume_oracle_kind.hpp"

namespace ttc
{

// Arithmetic of GenerateSamples and of the sample multiset X.
//   Auto (default):    None when the lattice coordinates k = x 10^b fit
//     comfortably in a double (b + log10 max|x| <= kAutoDoubleDigits), else
//     PointRepr.  The polytope coefficients are double from parsing on, so
//     within that range GMP adds no exactness (bit-identical estimates).
//   PointRepr (--gmp): the walk runs in double; lattice rounding, stored
//     points and the membership tests A x <= b are in GMP (mpf_float) at
//     `samplingPrecision` decimal digits.
//   None (--nogmp):    everything in double.
//   Full (--fullgmp):  the walk itself also runs in GMP.
enum class GmpMode
{
  Auto,
  PointRepr,
  None,
  Full
};

// Auto picks double while b + log10(max |x_i|) <= 15: the integer lattice
// coordinates are then exact (|k| <= 10^15 < 2^53 / 9).  Rounding x 10^b in
// double perturbs the fractional part that drives the randomized rounding by
// at most 10^15 2^-53 ~ 0.1 lattice spacing, half the resolution of the
// double walk's own proposals (ulp(x) 10^b), which GMP point storage shares.
constexpr double kAutoDoubleDigits = 15.0;

// b + log10 max(1, max|x|): the decimal digits of the largest lattice coordinate.
inline double latticeDigits(int b, double maxAbsCoordinate)
{
  return double(b) + std::log10(std::max(1.0, maxAbsCoordinate));
}

// GmpMode::Auto resolved for lattice coordinates of `digits` decimal digits.
inline GmpMode resolveGmpMode(GmpMode requested, double digits)
{
  if (requested != GmpMode::Auto) return requested;
  return digits <= kAutoDoubleDigits ? GmpMode::None : GmpMode::PointRepr;
}

// GMP decimal digits: b + 32, or the --precision minimum if larger.
inline int gmpDigits(int minimumDigits, int b) { return std::max(minimumDigits, b + 32); }

struct VolumeOptions
{
  double epsilon = 0.8;  // --epsilon, 0 < epsilon <= 1
  double delta = 0.2;    // --delta, 0 < delta < 1
  // Base seed of every random stream (seed_streams.hpp).  main.cpp sets it
  // from --seed only when the flag is given, so seedless runs keep 123.
  unsigned seed = 123;

  // VolEst (line 8): oracle, walk length (< 0: defaultOracleWalkLength) and
  // whether to round each polytope first (cleared by --no-round).
  VolumeOracle volumeOracle = kDefaultVolumeOracle;  // --vol-oracle
  int volumeWalkLength = -1;                         // --walklen-vol
  bool roundPolytopes = true;

  // GenerateSamples (line 19): walk and walk length (< 0: defaultSampleWalkLength).
  SamplerWalk samplerWalk = SamplerWalk::Billiard;  // --sampler
  int sampleWalkLength = -1;                        // --walklen-samp

  GmpMode gmpMode = GmpMode::Auto;  // --nogmp / --gmp / --fullgmp
  // Minimum GMP decimal digits, separate from the lattice precision b;
  // 0 means b + 32.
  int precision = 0;  // --precision

  bool cddSimplify = true;    // cleared by --no-cdd-simp
  std::string dumpInePrefix;  // --dump-ine PREFIX: write P_i to PREFIX_cubeN.ine
};

// Settings fixed once every polytope is canonicalized and has its balls,
// before the (expensive) volume oracles run.
struct VolumeSetup
{
  std::size_t dimension = 0;
  std::size_t keptPolytopes = 0;  // positive-volume polytopes after cdd
  unsigned volumeWalkLength = 0;
  unsigned sampleWalkLength = 0;
  int latticePrecision = 0;              // b: lattice spacing is 10^-b
  int samplingPrecision = 0;             // GMP decimal digits (0 when sampling in double)
  GmpMode gmpMode = GmpMode::PointRepr;  // resolved (never Auto)
  double latticeDigits = 0.0;            // b + log10(max |x_i| over all polytopes)
};

using SetupCallback = std::function<void(const VolumeSetup&)>;

}  // namespace ttc
