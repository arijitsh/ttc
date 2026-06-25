#pragma once

#include <cvc5/cvc5.h>
#include <utility>
#include <vector>

namespace arjun {

// Result of analyzing a projection set for redundant variables. The original
// projection set partitions into three disjoint groups:
//   * support  -- variables that must be kept (genuinely free given the rest)
//   * forced   -- variables pinned to a constant value in every model
//                 (backbone); paired with that value
//   * defined  -- variables implicitly defined by the support (Padoa /
//                 implicit definability), but not constant
//   * free     -- variables the formula does not constrain at all (they occur
//                 in no constraint), so both values are always admissible
// Removing forced/defined variables from the projection set preserves the
// *unweighted* projected count; removing a free variable halves it (each free
// variable doubles the model count). For weighted counting the caller must
// account for the weight contribution of each removed variable (see the
// parser).
struct Reduction
{
    std::vector<cvc5::Term> support;
    std::vector<std::pair<cvc5::Term, bool>> forced;
    std::vector<cvc5::Term> defined;
    std::vector<cvc5::Term> free;
};

// Classify the projection variables of `formula` into kept/forced/defined.
// The solver must be incremental; its assertion stack is restored on return.
Reduction analyzeProjectionSet(cvc5::Solver& solver,
                               const cvc5::Term& formula,
                               const std::vector<cvc5::Term>& projVars);

// Minimize the projection set by removing forced and implicitly definable
// variables. Updates projVars in place to contain an independent support of
// formula. (Thin wrapper over analyzeProjectionSet for callers that do not
// need the per-variable classification, e.g. unweighted counting.)
void minimizeProjectionSet(cvc5::Solver& solver,
                           const cvc5::Term& formula,
                           std::vector<cvc5::Term>& projVars);
}
