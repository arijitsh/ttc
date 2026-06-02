#ifndef D4_INTERFACE_HPP
#define D4_INTERFACE_HPP
#include <cstdint>
#include <cvc5/cvc5.h>
#include <vector>

// Count models of a CNF formula using the D4 projected model counter.
// `numVars` is the number of variables in the CNF. `cnf` contains one
// clause per entry, using positive integers for positive literals and
// negative integers for negated literals. `projVars` lists the variables
// to project on (1-based indices). `varOrder` provides a fixed branching
// order (lower values have higher priority) typically derived from a tree
// decomposition. `solver` is the original SMT solver that should follow
// D4's Boolean search, receiving assertions for each propagated literal.
// `idxToTerm` maps CNF variable indices back to their corresponding SMT
// terms.
struct D4Statistics {
  std::uint64_t cacheHits = 0;
  std::uint64_t cacheMisses = 0;
};

std::uint64_t d4Count(int numVars, const std::vector<std::vector<int>> &cnf,
                      const std::vector<int> &projVars,
                      const std::vector<double> &varOrder, cvc5::Solver &solver,
                      const std::vector<cvc5::Term> &idxToTerm, bool useSmtCache,
                      bool noDecompose, bool useResidualSimplifier,
                      bool fullDecompose, bool useMono, bool monoTrue,
                      bool useComponentCache, D4Statistics *stats);

#endif // D4_INTERFACE_HPP
