#include "pact/hash_constraint.hpp"

#include <sstream>
#include <utility>

#include "logger.hpp"

HashConstraint::HashConstraint(cvc5::Term fallback, std::vector<XorClause> clauses)
    : d_fallback(std::move(fallback)), d_xorClauses(std::move(clauses))
{
}

bool HashConstraint::isNull() const
{
  return d_fallback.isNull() && d_xorClauses.empty();
}

bool HashConstraint::hasXorClauses() const
{
  return !d_xorClauses.empty();
}

void HashConstraint::assertToSolver(cvc5::Solver& solver, bool useNativeXor) const
{
  bool canGoNative = useNativeXor && !d_xorClauses.empty();
  for (const XorClause& clause : d_xorClauses)
  {
    if (clause.terms.empty())
    {
      // An empty parity cannot be expressed as a native XOR clause; fall back
      // to the Boolean/bit-vector encoding for the whole hash.
      canGoNative = false;
      break;
    }
  }

  if (canGoNative)
  {
    const bool haveActivation = !d_activation.isNull();
    for (const XorClause& clause : d_xorClauses)
    {
      if (Trace.isEnabled("xor"))
      {
        Trace("xor") << "assertXorClause " << clauseToString(clause)
                      << std::endl;
      }
      // Hand the parity directly to the SAT solver (CaDiCaL) as a native XOR
      // clause so it is solved by Gauss-Jordan elimination rather than a
      // Tseitin CNF expansion. In --xor-activation literal we fold the hash's
      // indicator variable into every parity row: with it assigned false the
      // row enforces the real parity, and left free the row is vacuous, so the
      // hash can be toggled by an assumption without rebuilding the solver.
      if (haveActivation)
      {
        std::vector<cvc5::Term> terms = clause.terms;
        terms.push_back(d_activation);
        solver.assertXorClause(terms, clause.rhs);
      }
      else
      {
        solver.assertXorClause(clause.terms, clause.rhs);
      }
    }
    return;
  }

  if (!d_fallback.isNull())
  {
    solver.assertFormula(d_fallback);
  }
}

std::string HashConstraint::toString() const
{
  if (!d_xorClauses.empty())
  {
    std::ostringstream oss;
    for (std::size_t i = 0; i < d_xorClauses.size(); ++i)
    {
      if (i > 0)
      {
        oss << ", ";
      }
      oss << clauseToString(d_xorClauses[i]);
    }
    return oss.str();
  }

  std::ostringstream oss;
  oss << d_fallback;
  return oss.str();
}

std::string HashConstraint::clauseToString(const XorClause& clause) const
{
  std::ostringstream oss;
  oss << "(xor";
  for (const cvc5::Term& term : clause.terms)
  {
    oss << ' ' << term;
  }
  oss << ") = " << (clause.rhs ? "true" : "false");
  return oss.str();
}

