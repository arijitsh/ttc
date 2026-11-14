#pragma once

// Result rows of the LRA union-volume engine (kept free of cvc5/volesti so the
// union estimator and its tests can use them).

#include <cstddef>
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

struct VolumeComputationResult
{
  std::vector<VolumeComputationRow> rows;
  std::size_t finalSampleCount = 0;
  double volumeEstimate = 0.0;
  std::size_t totalSamplesGenerated = 0;
  std::size_t totalSamplesDeleted = 0;
  double volumeComputationTime = 0.0;
  // Inner ball (volesti IPM + Chebyshev LP) and outer-radius LPs, phase 1 and
  // the per-polytope sampler restarts; excluded from samplingTime.
  double ballTime = 0.0;
  double samplingTime = 0.0;
  // Union bookkeeping: membership of stored samples in each new polytope,
  // down-sampling and storing new samples (GMP unless --nogmp).
  double filterTime = 0.0;
  int samplingPrecision = 0;  // GMP decimal digits used for sample bookkeeping
  int latticePrecision = 0;   // b: lattice spacing is 10^-b
};

}  // namespace ttc
