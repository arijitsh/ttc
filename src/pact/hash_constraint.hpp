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

  // Optional activation ("indicator") variable for --xor-activation literal:
  // when set, it is folded into every native XOR parity so the whole hash can
  // be toggled by an assumption (assume it false -> the parity is enforced;
  // leave it free -> the parity row has a free GF(2) variable and is vacuous).
  // This lets a single solver hold the whole hash pool and switch hashes on/off
  // without the per-galloping-level solver rebuild the rebuild mode needs.
  void setActivation(cvc5::Term activation) { d_activation = std::move(activation); }
  const cvc5::Term& activation() const { return d_activation; }

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
  cvc5::Term d_activation;
};

