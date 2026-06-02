#ifndef VAR_ORDER_HPP
#define VAR_ORDER_HPP

#include <cvc5/cvc5.h>
#include <vector>

// Compute a heuristic ordering of projection variables based on
// a tree decomposition of the variable incidence graph.
// The returned vector contains the projection variables reordered
// so that variables with smaller treewidth appear earlier.
std::vector<cvc5::Term> computeProjVarOrder(const cvc5::Term& formula,
                                            const std::vector<cvc5::Term>& projVars,
                                            cvc5::Solver& solver,
                                            bool printTD = false,
                                            bool contract = true,
                                            bool netrel = false);

#endif // VAR_ORDER_HPP
