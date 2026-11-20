#pragma once

// The AllSAT engine behind ToDNF (Algorithm 1, line 2): the HALL circuit
// AllSAT solver of allsat-circuits (Fried-Nadel-Shalmon), run on the AIG F_bc.
// Each model it finds is generalized to a cube over the AIG inputs and then
// blocked, so the cubes together cover every model of the output but may
// overlap.  The solver's own types stay inside aig_allsat.cpp.

#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include "features.hpp"

namespace ttc
{

// One generalized model: (input variable, value) for the inputs in the order
// the solver reports them; an empty value is a don't-care.
using InputCube = std::vector<std::pair<std::uint32_t, std::optional<bool>>>;

std::vector<InputCube> allSatInputCubes(const BooleanAbstractionAigData& aig);

}  // namespace ttc
