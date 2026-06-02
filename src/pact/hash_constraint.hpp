#pragma once

#include <cvc5/cvc5.h>

#include <string>
#include <vector>

// Represents a single XOR clause consisting of a list of Boolean terms whose
// parity should equal the provided right-hand side.
struct XorClause
{
  std::vector<cvc5::Term> terms;
  bool rhs = false;
  cvc5::Term fallback;
};

// Encapsulates one randomly generated hash constraint. Each hash constraint may
// consist of one or more XOR clauses (for example, a Boolean clause and a
// bit-vector clause) together with a fallback term that encodes the full
// constraint when native XOR support is not available.
class HashConstraint
{
 public:
  HashConstraint() = default;
  HashConstraint(cvc5::Term fallback, std::vector<XorClause> clauses);

  bool isNull() const;
  bool hasXorClauses() const;

  const cvc5::Term& fallback() const { return d_fallback; }
  const std::vector<XorClause>& xorClauses() const { return d_xorClauses; }

  void assertToSolver(cvc5::Solver& solver, bool useNativeXor) const;
  std::string toString() const;

  bool operator==(const HashConstraint& other) const
  {
    return d_fallback == other.d_fallback;
  }

 private:
  std::string clauseToString(const XorClause& clause) const;

  cvc5::Term d_fallback;
  std::vector<XorClause> d_xorClauses;
};

