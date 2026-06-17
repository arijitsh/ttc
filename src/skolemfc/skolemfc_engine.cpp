#include "skolemfc/skolemfc_engine.hpp"

#include "features.hpp"
#include "logger.hpp"

#include <cstdlib>
#include <functional>
#include <unordered_set>

// SkolemFC library (links libskolemfc + CryptoMiniSat / ApproxMC / UniGen).
#include "skolemfc.h"

namespace ttc
{
namespace skolem
{

namespace
{

// Width in bits contributed by a projection term to the bit-blasted CNF.
std::uint32_t bitWidth(const cvc5::Term& term)
{
  cvc5::Sort sort = term.getSort();
  if (sort.isBitVector())
  {
    return sort.getBitVectorSize();
  }
  if (sort.isBoolean())
  {
    return 1;
  }
  return 0;
}

}  // namespace

SkolemFcResult countSkolemFunctions(cvc5::Solver& solver,
                                    const cvc5::Term& formula,
                                    const std::vector<cvc5::Term>& inputVars,
                                    const cvc5::Term& outputVar,
                                    std::uint64_t seed,
                                    double epsilon,
                                    double delta,
                                    int verbosity)
{
  SkolemFcResult result;

  // Bit-blast the formula. The sampling set is requested for the function
  // inputs first and the function output last, so getBitblastedCnf returns the
  // input bits as the leading, contiguous block of `samplingVars`.
  std::vector<cvc5::Term> projectionVars = inputVars;
  projectionVars.push_back(outputVar);

  ttc::BitblastedCnfData data;
  try
  {
    data = ttc::getBitblastedCnf(solver, formula, projectionVars);
  }
  catch (const std::exception& ex)
  {
    result.error = std::string("bit-blasting failed: ") + ex.what();
    return result;
  }

  // Number of leading sampling vars that belong to the function inputs.
  std::size_t inputBitCount = 0;
  for (const cvc5::Term& v : inputVars)
  {
    inputBitCount += bitWidth(v);
  }
  if (inputBitCount > data.samplingVars.size())
  {
    result.error = "internal error: fewer sampling bits than function inputs";
    return result;
  }

  // The forall set is the bits of the function inputs; the exists set is every
  // other CNF variable (function output bits and all internal bit-blast
  // variables), matching SkolemFC's QDIMACS convention (e = all non-a vars).
  std::unordered_set<std::uint32_t> forallSet;
  forallSet.reserve(inputBitCount * 2);
  for (std::size_t i = 0; i < inputBitCount; ++i)
  {
    forallSet.insert(data.samplingVars[i]);
  }

  std::size_t numClauses = 0;
  for (std::uint32_t v : data.clauses)
  {
    if (v == 0) ++numClauses;
  }

  // Exists variables (1-based): every CNF variable that is not a forall bit.
  std::vector<std::uint32_t> existsVars;
  existsVars.reserve(data.numVars - forallSet.size());
  for (std::uint32_t v = 1; v <= data.numVars; ++v)
  {
    if (forallSet.find(v) == forallSet.end()) existsVars.push_back(v);
  }

  result.numVars = data.numVars;
  result.numClauses = numClauses;
  result.numForall = forallSet.size();
  result.numExists = existsVars.size();

  Log(2) << "[skolemfc] bit-blasted CNF: " << data.numVars << " vars, "
         << numClauses << " clauses, " << forallSet.size() << " forall bits"
         << std::endl;

  // Helper: split the flat 0-separated DIMACS literal stream into clauses and
  // feed each to `emit` (with optional remapping of selected variables).
  auto forEachClause = [&](const std::function<CMSat::Lit(std::int32_t)>& mapLit,
                           const std::function<void(const std::vector<CMSat::Lit>&)>&
                               emit) {
    std::vector<CMSat::Lit> cl;
    for (std::uint32_t v : data.clauses)
    {
      if (v == 0)
      {
        emit(cl);
        cl.clear();
        continue;
      }
      cl.push_back(mapLit(static_cast<std::int32_t>(v)));
    }
    if (!cl.empty()) emit(cl);
  };
  auto plainLit = [](std::int32_t lit) {
    return CMSat::Lit(static_cast<std::uint32_t>(std::abs(lit)) - 1, lit < 0);
  };

  // (1) Is the relation F(input, output) satisfiable at all? If not there is no
  // realisable function and the Skolem-function count is 0.
  {
    CMSat::SATSolver sat;
    sat.set_verbosity(0);
    sat.new_vars(data.numVars);
    forEachClause(plainLit,
                  [&](const std::vector<CMSat::Lit>& cl) { sat.add_clause(cl); });
    if (sat.solve() == CMSat::l_False)
    {
      result.ok = true;
      result.satisfiable = false;
      result.log2Count = 0.0;
      Log(1) << "[skolemfc] relation is unsatisfiable; count is 0" << std::endl;
      return result;
    }
  }

  // (2) Determinism check. Build SkolemFC's "G formula" -- two copies of F that
  // share the forall (input) bits but use independent exists variables, plus a
  // gadget asserting the two exists assignments differ -- and SAT-check it. If
  // it is UNSAT every input has at most one output, so there is exactly one
  // total Skolem function (count 1). This mirrors create_g_formula() but solves
  // it with a plain SAT call rather than handing it to the (sampling) counter,
  // which aborts on an unsatisfiable instance.
  {
    // Primed copy of each exists variable, and one auxiliary var per exists.
    std::unordered_map<std::uint32_t, std::uint32_t> primeOf;
    primeOf.reserve(existsVars.size() * 2);
    std::uint32_t next = data.numVars;
    for (std::uint32_t e : existsVars) primeOf[e] = ++next;
    std::uint32_t firstAux = next;
    std::uint32_t totalVars = data.numVars + 2 * existsVars.size();

    CMSat::SATSolver g;
    g.set_verbosity(0);
    g.new_vars(totalVars);

    // F over the original variables.
    forEachClause(plainLit,
                  [&](const std::vector<CMSat::Lit>& cl) { g.add_clause(cl); });
    // F over the primed exists variables (forall bits are shared, untouched).
    auto primeLit = [&](std::int32_t lit) {
      std::uint32_t var = static_cast<std::uint32_t>(std::abs(lit));
      auto it = primeOf.find(var);
      std::uint32_t mapped = (it != primeOf.end()) ? it->second : var;
      return CMSat::Lit(mapped - 1, lit < 0);
    };
    forEachClause(primeLit,
                  [&](const std::vector<CMSat::Lit>& cl) { g.add_clause(cl); });
    // differ gadget: aux_e <-> (e != e'), and at least one aux_e true.
    std::vector<CMSat::Lit> anyDiff;
    anyDiff.reserve(existsVars.size());
    for (std::size_t i = 0; i < existsVars.size(); ++i)
    {
      std::uint32_t e = existsVars[i] - 1;            // 0-based
      std::uint32_t ep = primeOf[existsVars[i]] - 1;  // 0-based
      std::uint32_t aux = firstAux + i;               // 0-based (firstAux is 1-based count)
      // (e , e', ~aux) and (~e, ~e', ~aux): aux => e != e'.
      g.add_clause({CMSat::Lit(e, false), CMSat::Lit(ep, false),
                    CMSat::Lit(aux, true)});
      g.add_clause({CMSat::Lit(e, true), CMSat::Lit(ep, true),
                    CMSat::Lit(aux, true)});
      anyDiff.push_back(CMSat::Lit(aux, false));
    }
    if (!anyDiff.empty()) g.add_clause(anyDiff);

    if (existsVars.empty() || g.solve() == CMSat::l_False)
    {
      result.ok = true;
      result.deterministic = true;
      result.log2Count = 0.0;
      Log(1) << "[skolemfc] function is deterministic; count is 1" << std::endl;
      return result;
    }
  }

  // Genuinely non-deterministic: drive SkolemFC through its library interface.
  SkolemFC::SklFC skolemfc(epsilon, delta, static_cast<std::uint32_t>(seed),
                           static_cast<std::uint32_t>(verbosity));
  skolemfc.new_vars(data.numVars);
  forEachClause(plainLit, [&](const std::vector<CMSat::Lit>& cl) {
    skolemfc.add_clause(cl);
  });
  for (std::uint32_t v : forallSet) skolemfc.add_forall_var(v - 1);
  for (std::uint32_t v : existsVars) skolemfc.add_exists_var(v - 1);

  skolemfc.apply_default_config();
  skolemfc.count();

  if (skolemfc.count_succeeded())
  {
    result.ok = true;
    result.log2Count = skolemfc.get_count_log2();
  }
  else
  {
    result.error = "SkolemFC did not produce a count";
  }
  return result;
}

}  // namespace skolem
}  // namespace ttc
