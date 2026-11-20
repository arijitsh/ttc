#pragma once

// Line 2 of Algorithm 1:  F_bd <- ToDNF(F_bc).
//
// The AllSAT cubes of the AIG F_bc (aig_allsat.hpp), read through the map M
// as cubes over theory atoms.  Their polytopes P(c_i) cover the solutions of
// F and may overlap; Algorithm 1 measures their union.

#include <vector>

#include "features.hpp"
#include "volume/cube.hpp"

namespace ttc
{

// The cubes of F_bd, in the order AllSAT found them, after dropping
//   - cubes that give one atom both values (two inputs M maps to one atom),
//   - repeats of a cube already seen, and
//   - cubes whose literals strictly contain another cube's (their polytope
//     lies inside the other's, so they add nothing to the union).
std::vector<Cube> enumerateCubes(const BooleanAbstractionAigData& aig);

}  // namespace ttc
