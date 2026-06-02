#pragma once

#include <cvc5/cvc5.h>

#include <unordered_map>
#include <vector>

namespace d4 {

struct SimplifyResult {
  std::vector<cvc5::Term> assertions;
  bool unsat = false;
};

class ResidualSimplifier {
 public:
  explicit ResidualSimplifier(cvc5::Solver& solver);

  SimplifyResult simplify(const std::vector<cvc5::Term>& assertions);

 private:
  bool updateTrailSubstitutions();
  bool extractAssignment(const cvc5::Term& literal, cvc5::Term& term,
                         cvc5::Term& value) const;
  bool extractEqualityAssignment(const cvc5::Term& varCandidate,
                                 const cvc5::Term& valueCandidate,
                                 cvc5::Term& term, cvc5::Term& value) const;
  void flattenAnd(const cvc5::Term& t, std::vector<cvc5::Term>& out);

  cvc5::Solver& d_solver;
  std::vector<cvc5::Term> d_substTerms;
  std::vector<cvc5::Term> d_substValues;
  std::unordered_map<cvc5::Term, size_t> d_termToIndex;
  cvc5::Term d_trueTerm;
  cvc5::Term d_falseTerm;
};

}  // namespace d4
