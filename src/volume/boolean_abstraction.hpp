#pragma once

// Line 1 of Algorithm 1:  F_bc, M <- BoolAbs(F).
//
// F_bc is an And-Inverter Graph whose inputs are the theory atoms of F, and M
// (`BooleanAbstractionAigData::mapping`) sends each input variable to its
// atom.  ttc takes the AIG from its instrumented cvc5
// (getBooleanAbstractionAig in features.hpp) and falls back to
// buildBooleanAigFromTerm, which walks F itself, when that is unavailable.

#include <cvc5/cvc5.h>

#include "features.hpp"

namespace ttc
{

// AIG of `formula`: And, Or, Not, Implies, Ite, Xor and Boolean equality
// become gates; every other Boolean term is an atom (an AIG input).  Equal
// sub-terms and equal gates are shared, and constants are folded.
BooleanAbstractionAigData buildBooleanAigFromTerm(const cvc5::Term& formula);

}  // namespace ttc
