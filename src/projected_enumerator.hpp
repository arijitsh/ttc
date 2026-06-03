#pragma once

#include <cvc5/cvc5.h>

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "parser.hpp"

// Enumerates projected solutions by blocking each model on the projection
// variables and counting until UNSAT.
class ProjectedEnumerator
{
 public:
  ProjectedEnumerator(cvc5::Solver& solver,
                      const std::vector<cvc5::Term>& projectionVars);

  std::uint64_t count();
  long double countWeighted(
      const std::unordered_map<cvc5::Term, TTCParser::LiteralWeight>& weights);

  std::uint64_t getSmtCallCount() const { return d_smtCalls; }

 private:
  cvc5::Solver& d_solver;
  std::vector<cvc5::Term> d_projectionVars;
  std::uint64_t d_smtCalls = 0;
};
