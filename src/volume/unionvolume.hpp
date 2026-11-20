#pragma once

// The LRA volume engine: Algorithm 1 of the manuscript (kr25-ssm,
// sections/algo.tex), which estimates the volume of the solutions of an LRA
// formula F as the volume of a union of polytopes.  See unionvolume.cpp for
// the algorithm and README.md in this directory for the file of each step.

#include <vector>

#include <cvc5/cvc5.h>

#include "volume/cube.hpp"
#include "volume/volume_options.hpp"
#include "volume/volume_result.hpp"

namespace ttc
{

// Lines 3-22 of Algorithm 1 for the cubes F_bd of lines 1-2
// (TTCParser::computeCubes).  `realVariables` are the coordinates of R^n.
// onSetup is called once, before the volume oracles run, with the resolved
// settings; onRow once per kept polytope as the main loop finishes it.
// Throws std::invalid_argument on bad options and std::runtime_error when a
// step fails (unbounded polytope, every oracle run failing, ...).
VolumeComputationResult computeLraVolume(const std::vector<Cube>& cubes,
                                         const std::vector<cvc5::Term>& realVariables,
                                         const VolumeOptions& options = {},
                                         const RowCallback& onRow = {},
                                         const SetupCallback& onSetup = {});

}  // namespace ttc
