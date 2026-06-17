#ifndef TTC_SKOLEMFC_ENGINE_HPP
#define TTC_SKOLEMFC_ENGINE_HPP

#include <cvc5/cvc5.h>

#include <cstdint>
#include <string>
#include <vector>

namespace ttc
{
namespace skolem
{

// Result of counting the number of interpretations (Skolem functions) of an
// uninterpreted function that satisfy a formula.
struct SkolemFcResult
{
  bool ok = false;
  // log2 of the number of Skolem functions (the value SkolemFC reports on the
  // "s fc 2 ** <x>" line). The absolute count is 2 ** log2Count.
  double log2Count = 0.0;
  // The relation F(input, output) has at least one model. When false the
  // function cannot be realised at all, so the Skolem-function count is 0.
  bool satisfiable = true;
  // Every input admits at most one output (the "G formula" is UNSAT): there is
  // a single total Skolem function, so the count is exactly 1 (log2Count == 0).
  // Detected up front to avoid handing SkolemFC's sampler an unsatisfiable
  // formula, which it aborts on.
  bool deterministic = false;
  std::uint32_t numVars = 0;     // CNF variables after bit-blasting
  std::size_t numClauses = 0;    // CNF clauses
  std::size_t numForall = 0;     // forall (function input) bits
  std::size_t numExists = 0;     // exists (everything else, incl. output)
  std::string error;             // populated when ok == false
};

// Count the interpretations of an uninterpreted function by:
//   1. bit-blasting `formula` (the function application already replaced by
//      `outputVar`) into a model-preserving CNF via cvc5's getBitblastedCnf,
//   2. marking the bits of the function inputs (`inputVars`) as the universally
//      quantified (forall) variables and every remaining CNF variable -- the
//      function output bits and all internal bit-blast variables -- as the
//      existentially quantified (exists) variables,
//   3. handing the resulting QDIMACS instance to SkolemFC, whose count is the
//      number of valid input->output functions.
SkolemFcResult countSkolemFunctions(cvc5::Solver& solver,
                                    const cvc5::Term& formula,
                                    const std::vector<cvc5::Term>& inputVars,
                                    const cvc5::Term& outputVar,
                                    std::uint64_t seed,
                                    double epsilon,
                                    double delta,
                                    int verbosity);

}  // namespace skolem
}  // namespace ttc

#endif  // TTC_SKOLEMFC_ENGINE_HPP
