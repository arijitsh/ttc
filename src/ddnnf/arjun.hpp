#pragma once

#include <cvc5/cvc5.h>
#include <vector>

namespace arjun {
// Minimize the projection set by removing implicitly definable Boolean variables.
// The input solver must be incremental. The function updates projVars in place
// to contain an independent support of formula.
void minimizeProjectionSet(cvc5::Solver& solver,
                           const cvc5::Term& formula,
                           std::vector<cvc5::Term>& projVars);
}

