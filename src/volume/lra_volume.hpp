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
};

VolumeComputationResult computeLraVolume(
    const std::vector<Polytope>& polytopes,
    const std::vector<cvc5::Term>& realVariables,
    const std::function<void(const VolumeComputationRow&)>& onRow = {});

}  // namespace ttc

