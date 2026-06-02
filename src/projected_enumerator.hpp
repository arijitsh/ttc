#pragma once

#include <cvc5/cvc5.h>

#include <cstdint>
#include <vector>

// Enumerates projected solutions by blocking each model on the projection
// variables and counting until UNSAT.
class ProjectedEnumerator
{
 public:
  ProjectedEnumerator(cvc5::Solver& solver,
                      const std::vector<cvc5::Term>& projectionVars);

  std::uint64_t count();

  std::uint64_t getSmtCallCount() const { return d_smtCalls; }

 private:
  cvc5::Solver& d_solver;
  std::vector<cvc5::Term> d_projectionVars;
  std::uint64_t d_smtCalls = 0;
};
