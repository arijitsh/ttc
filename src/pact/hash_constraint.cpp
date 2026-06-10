#include "pact/hash_constraint.hpp"
#include <cstdint>
#include <cstdlib>

#include <sstream>
#include <utility>

#include "features.hpp"
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
  // --hash prime-gj: a mod-p row is enforced by cvc5's Z_p Gauss-Jordan
  // propagator, not bit-blasted. Hand the literals + weights + modulus straight
  // to assertModpClause; do NOT assert the (arithmetic) fallback term, so the
  // expensive multiply/urem is never built into CNF.
  for (const XorClause& clause : d_xorClauses)
  {
    if (clause.modulus != 0)
    {
      solver.assertModpClause(
          clause.terms, clause.weights, clause.rhsValue, clause.modulus);
      // The mod-p row lives in the CDCL(T) propagator, but the projection-bit
      // atoms (= (extract k k) #b1) otherwise appear in no asserted clause, so
      // the bit-vector theory never bit-blasts/links them to x's real bits and
      // the propagator reasons over free, disconnected copies. Force the theory
      // to bit-blast and link EACH bit-atom individually by asserting, per bit,
      // (= g_k lit_k) with a fresh unconstrained guard g_k. This bit-blasts and
      // ties lit_k to x's real bit (so the propagator's per-bit values match the
      // model) while adding no restriction on x -- and, unlike the arithmetic
      // fallback, blasts only the individual bits, never the multiply/mod. A
      // single parity link is insufficient: it would tie only the XOR of the
      // bits, letting individual atoms deviate and the mod-p sum drift.
      if (!std::getenv("TTC_MODP_NO_LINK"))
      {
        auto& tm = ttc::getTermBuilder(solver);
        static std::uint64_t s_guardCounter = 0;
        for (const cvc5::Term& lit : clause.terms)
        {
          cvc5::Term guard =
              tm.mkConst(tm.getBooleanSort(),
                         "__ttc_modp_link_" + std::to_string(s_guardCounter++));
          solver.assertFormula(tm.mkTerm(cvc5::Kind::EQUAL, {guard, lit}));
        }
      }
      return;
    }
  }

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

