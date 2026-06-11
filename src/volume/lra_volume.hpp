#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#include <cvc5/cvc5.h>

#include "volume/allsat_volume.hpp"

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
  double samplingTime = 0.0;
  int samplingPrecision = 0;  // GMP decimal digits used for sample bookkeeping
};

// Arithmetic backend for the sampling / union-volume bookkeeping.
//   PointRepr (default): the billiard walk runs in double, but sampled points
//     are stored and the membership tests (A x <= b) are evaluated in GMP
//     (boost mpf_float) at `samplingPrecision` decimal digits.
//   None (--nogmp):      the original all-double path (no GMP, no truncation).
//   Full (--fullgmp):    the billiard walk itself runs in GMP at the target
//     precision (patched volesti), then points are stored in GMP.
enum class GmpMode
{
  PointRepr,
  None,
  Full
};

// Tunable knobs for the LRA volume engine. A negative walk length means "use
// the volesti default" (volume: 10 + dim/10, sampling: 10).
struct VolumeOptions
{
  int volumeWalkLength = -1;  // --walklen-vol N
  int sampleWalkLength = -1;  // --walklen-samp N
  bool cddSimplify = true;    // cleared by --no-cdd-simp
  GmpMode gmpMode = GmpMode::PointRepr;
  // GMP precision (decimal digits) for sampling. 0 => auto, derived from the
  // polytopes via get_precision_from_cubes (src/cube_processor_nondis.py).
  int precision = 0;  // --precision N
};

VolumeComputationResult computeLraVolume(
    const std::vector<Polytope>& polytopes,
    const std::vector<cvc5::Term>& realVariables,
    const VolumeOptions& options = {},
    const std::function<void(const VolumeComputationRow&)>& onRow = {});

}  // namespace ttc

