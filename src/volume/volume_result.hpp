#pragma once

// Result of the LRA union-volume engine: one row per polytope of the main
// loop, the estimate, and where the time went.  Free of cvc5 and volesti so
// the union estimator and its tests can use it.

#include <cstddef>
#include <functional>
#include <vector>

namespace ttc
{

struct VolumeComputationRow
{
  std::size_t index = 0;
  double volume = 0.0;
  std::size_t samplesDeleted = 0;
  std::size_t totalSamples = 0;
  std::size_t samplesGenerated = 0;
};

using RowCallback = std::function<void(const VolumeComputationRow&)>;

struct VolumeComputationResult
{
  std::vector<VolumeComputationRow> rows;
  std::size_t finalSampleCount = 0;
  double volumeEstimate = 0.0;
  std::size_t totalSamplesGenerated = 0;
  std::size_t totalSamplesDeleted = 0;
  double volumeComputationTime = 0.0;
  // Inner balls and outer radii (LPs) of the polytopes and of the sampling
  // walks' bodies; excluded from samplingTime.
  double ballTime = 0.0;
  double samplingTime = 0.0;
  // Bookkeeping of X: membership tests, halvings and appends.
  double filterTime = 0.0;
  int samplingPrecision = 0;  // GMP decimal digits used for sample bookkeeping
  int latticePrecision = 0;   // b: lattice spacing is 10^-b
};

}  // namespace ttc
