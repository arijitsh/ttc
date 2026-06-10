#include "pact/pact.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <unordered_map>
#include <vector>
#include <cctype>

#include "features.hpp"
#include "profiler.hpp"

namespace
{
struct ProjectionSplit
{
  std::vector<cvc5::Term> all;
  std::vector<cvc5::Term> booleans;
  std::vector<cvc5::Term> bitvectors;
};

ProjectionSplit splitProjectionVars(const std::vector<cvc5::Term>& vars)
{
  ProjectionSplit split;
  split.all.reserve(vars.size());
  for (const cvc5::Term& t : vars)
  {
    if (t.getSort().isBoolean())
    {
      split.booleans.push_back(t);
      split.all.push_back(t);
    }
    else if (t.getSort().isBitVector())
    {
      split.bitvectors.push_back(t);
      split.all.push_back(t);
    }
  }
  return split;
}
}  // namespace

Pact::Pact(cvc5::Solver& solver,
           const std::vector<cvc5::Term>& projectionVars,
           std::uint64_t seed,
           bool useNativeXor,
           bool bvPact,
           double epsilon,
           double delta)
    : d_solver(solver),
      d_rng(static_cast<std::mt19937::result_type>(seed)),
      d_counter(solver),
      d_boolHash(solver),
      d_bvHash(solver),
      d_epsilon(epsilon),
      d_delta(delta)
{
  if (d_epsilon <= 0.0)
  {
    d_epsilon = 0.8;
  }
  if (d_delta <= 0.0 || d_delta >= 1.0)
  {
    d_delta = 0.2;
  }
  // ApproxMC's reuse_models is on by default; TTC_NO_REUSE turns it off so the
  // effect of the optimisation can be measured against the plain search.
  if (std::getenv("TTC_NO_REUSE"))
  {
    d_reuseModels = false;
  }
  if (const char* hm = std::getenv("TTC_HASH_MODE"))
  {
    d_useModpGj = (std::string(hm) == "prime-gj");
  }
  // --xor-activation {rebuild,literal}; main maps the flag to TTC_XOR_ACTIVATION.
  // Default 'rebuild': measured faster here than 'literal'. Although 'literal'
  // (one solver, hashes toggled by indicator assumptions) avoids the per-level
  // rebuild, the galloping search overshoots and every explored hash stays on
  // the solver -- on ttc's heavy count solver those accumulated, sampling-var-
  // coupled GF(2) rows cost more than the cheap rebuild they replace. (ApproxMC
  // tolerates it because its solver is a lean 45-var CNF, not a 183-var cvc5
  // re-encoding.) 'literal' is kept as an opt-in.
  if (const char* mode = std::getenv("TTC_XOR_ACTIVATION"))
  {
    d_xorActivationLiteral = (std::string(mode) == "literal");
  }
  else
  {
    d_xorActivationLiteral = false;
  }

  std::string logic;
  try
  {
    logic = d_solver.getLogic();
  }
  catch (const cvc5::CVC5ApiException&)
  {
    logic.clear();
  }
  if (!logic.empty())
  {
    std::string logicUpper = logic;
    std::transform(logicUpper.begin(),
                   logicUpper.end(),
                   logicUpper.begin(),
                   [](unsigned char ch) {
                     return static_cast<char>(std::toupper(ch));
                   });
    if (logicUpper != "ALL" && logicUpper.find("BV") == std::string::npos)
    {
      std::string newLogic = logic;
      std::size_t underscore = newLogic.find('_');
      if (underscore != std::string::npos)
      {
        newLogic.insert(underscore + 1, "BV");
      }
      else
      {
        newLogic = "BV" + newLogic;
      }
      try
      {
        d_solver.setLogic(newLogic);
      }
      catch (const cvc5::CVC5ApiException&)
      {
      }
    }
  }

  ProjectionSplit split = splitProjectionVars(projectionVars);
  d_projectionVars = std::move(split.all);
  d_booleanVars = std::move(split.booleans);
  d_bitvectorVars = std::move(split.bitvectors);

  d_booleanBvVars.clear();
  d_booleanBvVars.reserve(d_booleanVars.size());

  if (bvPact && !d_booleanVars.empty())
  {
    // BV-PACT mode: replace each Boolean projection variable with a fresh 1-bit
    // bit-vector variable and substitute (= bv #b1) for the Boolean throughout
    // the asserted formula.  The bit-vectors become the *primary* variables, so
    // the parity (XOR) hash constraints are genuine bit-vector constraints that
    // are routed to the bit-vector SAT backend (CryptoMiniSat) rather than being
    // folded back into the core SAT solver.
    auto& tm = ttc::getTermBuilder(d_solver);
    cvc5::Sort bvSort = tm.mkBitVectorSort(1);
    cvc5::Term bvOne = tm.mkBitVector(1, 1);

    std::vector<cvc5::Term> from;
    std::vector<cvc5::Term> to;
    std::vector<cvc5::Term> eqTerms;
    from.reserve(d_booleanVars.size());
    to.reserve(d_booleanVars.size());
    eqTerms.reserve(d_booleanVars.size());
    for (std::size_t i = 0; i < d_booleanVars.size(); ++i)
    {
      std::string name = "__ttc_pbv_" + std::to_string(i);
      cvc5::Term bvVar = tm.mkConst(bvSort, name);
      cvc5::Term eq = tm.mkTerm(cvc5::Kind::EQUAL, {bvVar, bvOne});
      from.push_back(d_booleanVars[i]);
      to.push_back(eq);
      d_booleanBvVars.push_back(bvVar);
      eqTerms.push_back(eq);
    }

    // Re-assert the formula with the Boolean projection variables substituted
    // out.  resetAssertions() keeps options and declarations intact.
    std::vector<cvc5::Term> assertions = d_solver.getAssertions();
    d_solver.resetAssertions();
    for (cvc5::Term& assertion : assertions)
    {
      d_solver.assertFormula(assertion.substitute(from, to));
    }

    // Count and hash over the fresh bit-vector projection variables.
    d_projectionVars.clear();
    d_projectionVars.reserve(d_booleanBvVars.size() + d_bitvectorVars.size());
    d_projectionVars.insert(
        d_projectionVars.end(), d_booleanBvVars.begin(), d_booleanBvVars.end());
    d_projectionVars.insert(
        d_projectionVars.end(), d_bitvectorVars.begin(), d_bitvectorVars.end());

    // BoolHash builds a ~0.5-density parity over the bit-vector parity terms; we
    // pass the (= bv #b1) literals for the (unused) native-XOR representation and
    // the bit-vectors themselves for the asserted bit-vector parity constraint.
    d_boolHash.setVariables(eqTerms, d_booleanBvVars);
    d_bvHash.setVariables(d_bitvectorVars);
  }
  else
  {
    // Boolean-PACT (backup) mode: keep the Boolean projection variables primary
    // and slave a 1-bit bit-vector to each one via (= bv (ite bool #b1 #b0)).
    if (!d_booleanVars.empty())
    {
      auto& tm = ttc::getTermBuilder(d_solver);
      cvc5::Sort bvSort = tm.mkBitVectorSort(1);
      cvc5::Term bvOne = tm.mkBitVector(1, 1);
      cvc5::Term bvZero = tm.mkBitVector(1, 0);
      for (std::size_t i = 0; i < d_booleanVars.size(); ++i)
      {
        std::string name = "__ttc_bool_proj_bv_" + std::to_string(i);
        cvc5::Term bvVar = tm.mkConst(bvSort, name);
        cvc5::Term ite = tm.mkTerm(cvc5::Kind::ITE, {d_booleanVars[i], bvOne, bvZero});
        cvc5::Term eq = tm.mkTerm(cvc5::Kind::EQUAL, {bvVar, ite});
        d_solver.assertFormula(eq);
        d_booleanBvVars.push_back(bvVar);
      }
    }
    d_boolHash.setVariables(d_booleanVars, d_booleanBvVars);
    d_bvHash.setVariables(d_bitvectorVars);
  }

  d_counter.setProjectionVars(&d_projectionVars);
  d_counter.setUseNativeXor(useNativeXor);
  d_useNativeXor = useNativeXor;

  // The galloping search needs an upper sentinel on the number of useful hash
  // constraints. Each parity hash roughly halves the model space, so no more
  // than the total number of projection bits can ever help: Booleans contribute
  // one bit each and bit-vectors their width.
  d_maxHashes = 0;
  for (const cvc5::Term& var : d_projectionVars)
  {
    cvc5::Sort sort = var.getSort();
    if (sort.isBoolean())
    {
      d_maxHashes += 1;
    }
    else if (sort.isBitVector())
    {
      d_maxHashes += sort.getBitVectorSize();
    }
  }
  if (d_maxHashes == 0)
  {
    d_maxHashes = 1;
  }

  // Word-level bit-vector hash families (--hash prime/lemire) split the space
  // into more than two cells per hash. They only make sense over a purely
  // bit-vector projection (main rejects mixing them with Boolean projection
  // variables); when active, fewer hashes are needed and each measurement scales
  // by the family's per-hash multiplier rather than 2. Otherwise the default XOR
  // multiplier of 2 is kept.
  if (!d_bitvectorVars.empty() && d_booleanVars.empty())
  {
    d_hashMultiplier = d_bvHash.perHashMultiplier();
  }
  if (d_hashMultiplier > 2.0)
  {
    const double bitsPerHash = std::log2(d_hashMultiplier);
    if (bitsPerHash > 0.0)
    {
      std::size_t refined = static_cast<std::size_t>(
                                std::ceil(d_maxHashes / bitsPerHash))
                            + 1;
      d_maxHashes = std::max<std::size_t>(refined, 1);
    }
  }

  d_baseAssertions = d_solver.getAssertions();
}

void Pact::rebuildCountSolver()
{
  // Build a fresh solver on the same term manager so the previously asserted
  // hash terms (and their accumulated bit-blasting in the SAT backend) are
  // discarded. Reusing one solver makes each round progressively slower.
  cvc5::TermManager& tm = d_solver.getTermManager();
  d_countSolver.emplace(tm);
  cvc5::Solver& cs = *d_countSolver;

  if (d_useNativeXor)
  {
    // Must be set before the solver is otherwise initialised. The native XOR is
    // asserted over the SAT literals of the projection booleans. cvc5
    // preprocessing substitutes top-level equality definitions (e.g.
    // '(= b (>= x 0))' -> b := (>= x 0)), which removes b's own literal from the
    // formula; the XOR then constrains b's now-free literal while the reported
    // model reconstructs b from its definition. They desync and the XOR is
    // silently not enforced on the counted assignment. Disabling this
    // substitution keeps the projection booleans as real variables.
    if (const char* spec = std::getenv("TTC_PP"))
    {
      // TTC_PP="opt1=val1,opt2=val2" for experimenting with preprocessing opts.
      std::string s(spec), cur;
      std::vector<std::string> parts;
      for (char c : s) { if (c == ',') { parts.push_back(cur); cur.clear(); } else cur += c; }
      if (!cur.empty()) parts.push_back(cur);
      for (const std::string& p : parts)
      {
        auto eq = p.find('=');
        if (eq == std::string::npos) continue;
        try { cs.setOption(p.substr(0, eq), p.substr(eq + 1)); }
        catch (const cvc5::CVC5ApiException&) {}
      }
    }
    else if (!std::getenv("TTC_NO_SIMP_NONE"))
    {
      try { cs.setOption("simplification", "none"); }
      catch (const cvc5::CVC5ApiException&) {}
    }
    try { cs.setOption("simplification-bcp", "false"); }
    catch (const cvc5::CVC5ApiException&) {}
  }

  // Mirror the solving-relevant options from the original solver so the count
  // solver behaves identically (logic, incrementality, model production and the
  // bit-vector / XOR backend selection set up by main).
  static constexpr const char* kOptions[] = {
      "print-success",  "incremental",        "produce-models",
      "bv-sat-solver",  "bv-to-bool",         "sat-use-native-xor",
      "sat-solver",     "portfolio-jobs",
  };
  for (const char* name : kOptions)
  {
    try
    {
      cs.setOption(name, d_solver.getOption(name));
    }
    catch (const cvc5::CVC5ApiException&)
    {
    }
  }

  if (!d_useNativeXor && std::getenv("TTC_XOR_CNF"))
  {
    // --xorcnf baseline: hash with the Boolean parity *formula* (the fallback
    // path, not assertXorClause) and force CaDiCaL with native XOR disabled, so
    // cvc5's CNF stream Tseitin-expands each XOR into ordinary clauses. The
    // parity is then solved as plain CNF by CaDiCaL -- no Gauss-Jordan at all --
    // which is the natural baseline for the XOR-reasoning comparison.
    try
    {
      cs.setOption("sat-solver", "cadical");
      cs.setOption("sat-use-native-xor", "false");
    }
    catch (const cvc5::CVC5ApiException&)
    {
    }
  }

  if (d_useNativeXor)
  {
    // Native XOR clauses (assertXorClause) are only supported when the core SAT
    // solver is CaDiCaL; Minisat aborts on them. Force CaDiCaL on the counting
    // solver so the parity hashes are solved by Gauss-Jordan elimination.
    try
    {
      cs.setOption("sat-solver", "cadical");
      cs.setOption("sat-use-native-xor", "true");
    }
    catch (const cvc5::CVC5ApiException&)
    {
    }
    try
    {
      cs.setXorAssertionVerbose(false);
    }
    catch (const cvc5::CVC5ApiException&)
    {
    }
  }

  // --hash prime-gj: the mod-p Gauss-Jordan rows live in the CaDiCaL CDCL(T)
  // propagator. The projection bits are tied to x via the per-bit linkage
  // clauses asserted in HashConstraint::assertToSolver; here we force the count
  // solver to CaDiCaL (for the propagator) and disable simplification so those
  // linkage clauses survive.
  if (const char* hm = std::getenv("TTC_HASH_MODE"))
  {
    if (std::string(hm) == "prime-gj")
    {
      // The per-bit linkage (= g_k lit_k) asserted alongside each mod-p row has
      // unconstrained guards; cvc5's nonclausal simplification would substitute
      // them away and drop the clauses, un-linking the bits from x. Disable
      // simplification so the linkage survives.
      try
      {
        cs.setOption("simplification", "none");
      }
      catch (const cvc5::CVC5ApiException&)
      {
      }
      // The mod-p Gauss-Jordan engine lives in the CaDiCaL CDCL(T) propagator,
      // so the count solver must be CaDiCaL (no native XOR -- prime-gj adds no
      // XOR clauses, and native-XOR-on-BV is unstable).
      try
      {
        cs.setOption("sat-solver", "cadical");
        cs.setOption("sat-use-native-xor", "false");
      }
      catch (const cvc5::CVC5ApiException&)
      {
      }
    }
  }
  try
  {
    std::string logic = d_solver.getLogic();
    if (!logic.empty())
    {
      cs.setLogic(logic);
    }
  }
  catch (const cvc5::CVC5ApiException&)
  {
  }

  for (const cvc5::Term& assertion : d_baseAssertions)
  {
    cs.assertFormula(assertion);
  }

  d_counter.setSolver(cs);
  d_counter.resetCache();
  // The fresh solver carries none of the pool's hashes yet (literal activation).
  d_assertedHashes = 0;
}

double Pact::calcErrorBound(std::size_t t, double p) const
{
  // sum_{k=ceil(t/2)}^{t} C(t,k) p^k (1-p)^(t-k), computed by walking down from
  // the k=t term so the binomial coefficients update multiplicatively (no
  // factorials / overflow). Verbatim from ApproxMC's Counter::calc_error_bound.
  double curr = std::pow(p, static_cast<double>(t));
  double sum = curr;
  const std::int64_t lo = static_cast<std::int64_t>(
      std::ceil(static_cast<double>(t) / 2.0));
  for (std::int64_t k = static_cast<std::int64_t>(t) - 1; k >= lo; --k)
  {
    curr *= static_cast<double>(k + 1) / static_cast<double>(t - k)
            * (1.0 - p) / p;
    sum += curr;
  }
  return sum;
}

Pact::Parameters Pact::getParameters() const
{
  constexpr double appmc7EpsCutoff = 5.0;
  const bool useAppmc7 = d_epsilon >= appmc7EpsCutoff;
  double alpha = -1.0;
  double beta = -1.0;

  // ApproxMC6 threshold (CAV23, Counter::set_up_probs_threshold_measurements,
  // dense / non-sparse so thresh_factor = 1).
  std::size_t pivot = static_cast<std::size_t>(
      1.0
      + 9.84 * (1.0 + 1.0 / d_epsilon) * (1.0 + 1.0 / d_epsilon)
            * (1.0 + d_epsilon / (1.0 + d_epsilon)));
  if (useAppmc7)
  {
    // ApproxMC7 Algorithm 4, line 6: for large epsilon the threshold is zero
    // and each measurement searches for the last non-empty cell.
    pivot = 0;
    beta = (1.0
            + std::sqrt(1.0 + 2.0 * (1.0 + d_epsilon) * (1.0 + d_epsilon)))
           / 2.0;
    alpha = beta - 1.0;
  }
  if (pivot == 0)
  {
    if (!useAppmc7)
    {
      pivot = 1;
    }
  }
  if (const char* e = std::getenv("TTC_PIVOT"))
  {
    pivot = static_cast<std::size_t>(std::atoll(e));
  }

  // Upper bounds on the per-measurement under-/over-estimation probabilities
  // (ApproxMC6, Lemma 4 of the CAV23 rounding paper). The epsilon breakpoints
  // pair with the rounding applied in roundCount().
  double pL;
  if (d_epsilon < std::sqrt(2.0) - 1.0)
  {
    pL = 0.262;
  }
  else if (d_epsilon < 1.0)
  {
    pL = 0.157;
  }
  else if (d_epsilon < 3.0)
  {
    pL = 0.085;
  }
  else if (d_epsilon < 4.0 * std::sqrt(2.0) - 1.0)
  {
    pL = 0.055;
  }
  else
  {
    pL = 0.023;
  }
  double pU = (d_epsilon < 3.0) ? 0.169 : 0.044;
  if (useAppmc7)
  {
    // ApproxMC7 Lemma 5.
    pL = 1.0 / (1.0 + alpha);
    pU = 1.0 / beta;
  }

  // Smallest odd number of measurements whose combined error bound <= delta
  // (median amplification, ApproxMC6 Algorithm 6).
  std::size_t iterations = 1;
  while (calcErrorBound(iterations, pL) + calcErrorBound(iterations, pU)
         > d_delta)
  {
    iterations += 2;
  }
  if (const char* it = std::getenv("TTC_ITERATIONS"))
  {
    iterations = static_cast<std::size_t>(std::atoll(it));
    if (iterations == 0)
    {
      iterations = 1;
    }
  }
  return {pivot, iterations, useAppmc7, alpha, beta};
}

double Pact::roundCount(double cellCount) const
{
  // ApproxMC6 rounding (CAV23 Algorithm 5): round the cell count up to an
  // epsilon-dependent multiple of the pivot. The pivot here is the un-rounded
  // 9.84*(1+1/eps)^2 (NOT the +1, *thresh_factor threshold above).
  const double pivot =
      9.84 * (1.0 + 1.0 / d_epsilon) * (1.0 + 1.0 / d_epsilon);
  if (d_epsilon < std::sqrt(2.0) - 1.0)
  {
    const double floor = std::sqrt(1.0 + 2.0 * d_epsilon) / 2.0 * pivot;
    return std::max(cellCount, floor);
  }
  if (d_epsilon < 1.0)
  {
    return std::max(cellCount, pivot / std::sqrt(2.0));
  }
  if (d_epsilon < 3.0)
  {
    return std::max(cellCount, pivot);
  }
  if (d_epsilon < 4.0 * std::sqrt(2.0) - 1.0)
  {
    return pivot;
  }
  return std::sqrt(2.0) * pivot;
}

HashConstraint Pact::generateHashConstraint()
{
  auto& tm = ttc::getTermBuilder(d_solver);
  std::vector<XorClause> clauses;
  clauses.reserve(2);

  if (auto boolClause = d_boolHash.randomClause(d_rng))
  {
    clauses.push_back(std::move(*boolClause));
  }
  if (auto bvClause = d_bvHash.randomClause(d_rng))
  {
    clauses.push_back(std::move(*bvClause));
  }

  cvc5::Term fallback;
  if (clauses.empty())
  {
    fallback = tm.mkBoolean(true);
  }
  else if (clauses.size() == 1)
  {
    fallback = clauses.front().fallback;
  }
  else
  {
    std::vector<cvc5::Term> parts;
    parts.reserve(clauses.size());
    for (const XorClause& clause : clauses)
    {
      parts.push_back(clause.fallback);
    }
    fallback = tm.mkTerm(cvc5::Kind::AND, parts);
  }

  HashConstraint hash(fallback, std::move(clauses));

  // In --xor-activation literal, give the hash a fresh indicator variable folded
  // into its parity (see HashConstraint::setActivation). Only useful for the
  // native-XOR path; the fallback Boolean encoding still rebuilds.
  if (d_xorActivationLiteral && d_useNativeXor && hash.hasXorClauses())
  {
    cvc5::Term act =
        tm.mkConst(tm.getBooleanSort(),
                   "xact_" + std::to_string(d_activationCounter++));
    hash.setActivation(act);
  }

  return hash;
}

Pact::ProbeResult Pact::probeHashCount(std::size_t hashCount)
{
  Parameters params = getParameters();
  const std::size_t threshold = params.appmc7 ? 1 : params.threshold;
  d_counter.resetStatistics();
  d_savedModels.clear();
  rebuildCountSolver();

  std::vector<HashConstraint> hashPool;
  hashPool.reserve(hashCount);
  while (hashPool.size() < hashCount)
  {
    hashPool.push_back(generateHashConstraint());
  }

  std::optional<std::size_t> attempt = d_counter.count(hashPool, threshold);
  if (d_progressCallback)
  {
    d_progressCallback({hashCount, attempt, "probe", 0, 0, threshold});
  }
  ProbeResult result;
  result.hashCount = hashCount;
  result.saturated = !attempt.has_value();
  result.cellCount = attempt.value_or(threshold);
  d_lastWinningHash = hashCount;
  d_lastWinningCellCount = result.cellCount;
  return result;
}

double Pact::estimateFromMeasurement(std::size_t hashCount,
                                     std::size_t cellCount) const
{
  Parameters params = getParameters();
  double cell = static_cast<double>(cellCount);
  if (params.appmc7 && hashCount > 0)
  {
    cell *= std::sqrt(2.0 * params.alpha / params.beta);
  }
  else if (hashCount > 0 && d_hashMultiplier == 2.0)
  {
    // ApproxMC6 rounding is calibrated for factor-2 hashing, where the winning
    // cell count always lands within [pivot/2, pivot] so the rounding floor is a
    // sound correction. The word-level families (--hash prime/lemire) divide the
    // space by their (much larger) multiplier per hash, so the winning cell
    // count can be as small as pivot/multiplier; applying the factor-2 floor
    // there grossly overestimates. For them use the raw cell count, i.e. the
    // classic SMTApproxMC estimator cellCount * multiplier^hashCount.
    cell = roundCount(cell);
  }
  // cell * d_hashMultiplier^hashCount. For the default XOR family
  // (d_hashMultiplier == 2) this is exactly std::ldexp(cell, hashCount); the
  // word-level families use their larger per-hash multiplier instead.
  return cell * std::pow(d_hashMultiplier, static_cast<double>(hashCount));
}

std::uint64_t Pact::count()
{
  Parameters params = getParameters();
  if (d_iterationOverride.has_value())
  {
    params.iterations = std::max<std::size_t>(1, *d_iterationOverride);
  }
  d_counter.resetStatistics();
  Trace("pact") << "[pact] Starting approximate count with pivot "
                 << params.threshold << " and " << params.iterations
                 << " iterations" << std::endl;

  double startTime = Log.elapsed();
  double lastProgress = -progress_interval;
  std::size_t currentRound = 0;
  // Peak number of parity hashes added at any evaluation across the whole
  // galloping search (reported as 'c max hashes:' for comparison with other
  // counters' XOR depth).
  std::size_t maxHashesUsed = 0;

  // Emits one row for every saturating-count evaluation (one per hash count
  // tried by the galloping search). `nextHash` is where the search will jump
  // next.
  auto report = [&](std::size_t hashCount,
                    const std::optional<std::size_t>& sols,
                    std::int64_t nextHash,
                    std::size_t reuseSat,
                    std::size_t reuseChecked) {
    if (hashCount > maxHashesUsed)
    {
      maxHashesUsed = hashCount;
    }
    if (d_progressCallback)
    {
      d_progressCallback(
          {hashCount,
           sols,
           std::to_string(nextHash),
           reuseSat,
           reuseChecked,
           params.appmc7 ? 1 : params.threshold});
    }
    if (Log.getVerbosity() == 0)
    {
      return;
    }
    if (d_quiet)
    {
      return;
    }
    double elapsed = Log.elapsed() - startTime;
    lastProgress = elapsed;
    std::ostringstream countStream;
    if (sols.has_value())
    {
      countStream << *sols;
    }
    else
    {
      countStream << ">=" << (params.appmc7 ? 1 : params.threshold);
    }
    // reuse_models tally: of `reuseChecked` stored solutions examined this
    // evaluation, `reuseSat` still satisfied the active hashes and were reused
    // without an SMT call. '-' when nothing was examined.
    std::ostringstream reuseStream;
    if (reuseChecked == 0)
    {
      reuseStream << '-';
    }
    else
    {
      reuseStream << reuseSat << '/' << reuseChecked;
    }
    std::ostringstream line;
    line << std::fixed << std::setprecision(2);
    line << "c " << std::setw(7) << elapsed;
    line << ' ' << std::setw(3) << currentRound;
    line << ' ' << std::setw(4) << hashCount;
    line << ' ' << std::setw(8) << countStream.str();
    line << ' ' << std::setw(9) << reuseStream.str();
    line << ' ' << std::setw(5) << nextHash;
    std::cout << line.str() << std::endl;
  };

  if (!params.appmc7 && d_assumeBaseSaturated)
  {
    d_baseSaturated = true;
  }
  else if (!params.appmc7)
  {
    // Base count with no hashing. If the formula has fewer than `threshold`
    // models in total the count is exact and we are done -- this mirrors
    // ApproxMC "counting without XORs". ApproxMC7's threshold is zero and its
    // level-0 check is part of appmc7_one_measurement_count instead.
    rebuildCountSolver();
    std::vector<HashConstraint> empty;
    std::optional<std::size_t> base = d_counter.count(empty, params.threshold);
    report(0, base, 0, 0, 0);
    if (base.has_value())
    {
      Trace("pact") << "[pact] Base model count succeeded without hashing: "
          << *base << std::endl;
      d_lastWinningHash = 0;
      d_lastWinningCellCount = *base;
      if (!d_quiet)
      {
        std::cout << "c max hashes: " << maxHashesUsed << std::endl;
      }
      return static_cast<std::uint64_t>(*base);
    }
    Trace("pact")
        << "[pact] Base count reached saturation; introducing random hashes"
        << std::endl;
    // Level 0 (no hashes) is now known to saturate; oneMeasurement reuses this
    // instead of re-enumerating the hash-free formula every measurement.
    d_baseSaturated = true;
  }
  else
  {
    d_baseSaturated = false;
  }

  std::vector<double> estimates;
  estimates.reserve(params.iterations);

  // Carries the winning hash count between measurements (ApproxMC's
  // prev_measure) so later rounds start their galloping search near the answer.
  std::int64_t prevMeasure = d_initialPrevMeasure;

  for (std::size_t iter = 0; iter < params.iterations; ++iter)
  {
    currentRound = iter + 1;
    Trace("pact") << "[pact] Iteration " << (iter + 1) << " of "
        << params.iterations << std::endl;

    // Rebuild the counting solver from scratch so each round solves at a stable
    // speed instead of degrading as hash clauses accumulate over time.
    rebuildCountSolver();

    // Each measurement draws a fresh pool of random parity hashes; within the
    // measurement the pool is reused so re-evaluating the same hash count is
    // deterministic. The saved-model store is tied to that pool, so it is
    // cleared here (ApproxMC's hm.clear() between measurements).
    std::vector<HashConstraint> hashPool;
    d_savedModels.clear();
    MeasurementResult m =
        params.appmc7
            ? oneMeasurementAppmc7(iter, prevMeasure, hashPool, report)
            : oneMeasurement(
                  iter, prevMeasure, hashPool, params.threshold, report);

    // ApproxMC6 rounds the cell count; ApproxMC7 applies its alpha/beta
    // adjustment. Both happen before scaling by 2^hashCount. With no hashes the
    // cell count is exact and must not be inflated.
    double approx = estimateFromMeasurement(m.hashCount, m.cellCount);
    d_lastWinningHash = m.hashCount;
    d_lastWinningCellCount = m.cellCount;
    Trace("pact") << "[pact] Iteration " << (iter + 1) << " estimate: "
        << approx << " (" << m.cellCount << " @ " << m.hashCount << ")"
        << std::endl;
    estimates.push_back(approx);

    if (Log.getVerbosity() != 0 && !d_quiet)
    {
      double elapsed = Log.elapsed() - startTime;
      std::ostringstream line;
      line << std::fixed << std::setprecision(2);
      line << "c " << std::setw(7) << elapsed;
      line << ' ' << std::setw(3) << currentRound;
      line << ' ' << std::setw(4) << m.hashCount;
      line << ' ' << std::setw(8) << m.cellCount;
      line << ' ' << std::setw(9) << "-";
      line << ' ' << std::setw(5) << "win";
      line << ' ' << std::setw(13)
           << static_cast<std::uint64_t>(std::llround(approx));
      std::cout << line.str() << std::endl;
      lastProgress = elapsed;
    }
  }

  std::sort(estimates.begin(), estimates.end());
  double medianEstimate =
      estimates.empty() ? 0.0 : estimates[estimates.size() / 2];
  std::uint64_t result = static_cast<std::uint64_t>(std::llround(medianEstimate));
  Trace("pact") << "[pact] Median estimate across iterations: " << result
      << std::endl;
  if (!d_quiet)
  {
    std::cout << "c max hashes: " << maxHashesUsed << std::endl;
  }
  if (std::getenv("TTC_PACT_STATS") != nullptr)
  {
    try
    {
      cvc5::Solver& statSolver = d_countSolver ? *d_countSolver : d_solver;
      std::string stats = statSolver.getStatistics().toString();
      std::istringstream iss(stats);
      std::string line;
      while (std::getline(iss, line))
      {
        if (line.find("cryptominisat") != std::string::npos
            || line.find("BVSolverBitblast") != std::string::npos
            || line.find("cadical") != std::string::npos
            || line.find("gauss") != std::string::npos
            || line.find("xor") != std::string::npos
            || line.find("bitblast") != std::string::npos)
        {
          std::cerr << "[stat] " << line << std::endl;
        }
      }
    }
    catch (const std::exception& ex)
    {
      std::cerr << "[stat] unavailable: " << ex.what() << std::endl;
    }
  }
  return result;
}

Pact::MeasurementResult Pact::oneMeasurement(
    std::size_t iter,
    std::int64_t& prevMeasure,
    std::vector<HashConstraint>& hashPool,
    std::size_t threshold,
    const ReportFn& report)
{
  // Galloping search over the number of hash constraints, mirroring ApproxMC's
  // one_measurement_count. The two sentinels lowerFib (loIndex) and upperFib
  // (hiIndex) always bracket the answer: lowerFib is the largest hash count
  // known to saturate the pivot, upperFib the smallest known to fall below it.
  // We search exponentially until the bracket is small enough, then binary.
  const std::int64_t totalMaxHashes = static_cast<std::int64_t>(d_maxHashes);
  std::int64_t lowerFib = 0;
  std::int64_t upperFib = totalMaxHashes;
  std::int64_t numExplored = 0;
  std::int64_t hashCnt = prevMeasure;
  std::int64_t hashPrev = hashCnt;

  // thresholdSols[h] == true  : `h` hashes saturated the pivot (>= threshold)
  // thresholdSols[h] == false : `h` hashes gave an exact count (< threshold)
  // solsForHash[h]            : the exact count recorded when below the pivot
  std::unordered_map<std::int64_t, bool> thresholdSols;
  std::unordered_map<std::int64_t, std::int64_t> solsForHash;
  thresholdSols[totalMaxHashes] = false;
  solsForHash[totalMaxHashes] = 1;

  auto activeHashes = [&](std::int64_t level) {
    while (static_cast<std::int64_t>(hashPool.size()) < level)
    {
      hashPool.push_back(generateHashConstraint());
    }
    return std::vector<HashConstraint>(hashPool.begin(),
                                       hashPool.begin() + level);
  };

  while (numExplored < totalMaxHashes)
  {
    if (hashCnt < 0)
    {
      hashCnt = 0;
    }
    const std::int64_t curHashCnt = hashCnt;
    std::optional<std::size_t> attempt;
    // reuse_models tally for this evaluation, threaded into report().
    std::size_t curReuseSat = 0;
    std::size_t curReuseChecked = 0;
    if (curHashCnt == 0 && d_baseSaturated)
    {
      // Level 0 is the hash-free formula, already known to saturate the pivot
      // (see d_baseSaturated). Skip the rebuild + full re-enumeration: on large
      // formulas with a small count the search returns to level 0 every
      // measurement, and that base enumeration is the dominant cost.
      Trace("pact") << "[pact] Reusing saturated base count at level 0"
          << std::endl;
    }
    else if (d_xorActivationLiteral && d_useNativeXor)
    {
      // Literal activation: the whole hash pool lives on one solver. Generate
      // hashes up to hashCnt, assert any not yet present (once each, with their
      // indicator folded into the parity), then count with the first hashCnt
      // indicators assumed false -- so exactly those hashes are enforced and the
      // rest stay vacuous. No per-level rebuild; the search moving down just
      // assumes fewer indicators.
      activeHashes(hashCnt);  // grow hashPool to >= hashCnt
      cvc5::Solver& cs = d_countSolver ? *d_countSolver : d_solver;
      Trace("pact") << "[pact] Evaluating hash count " << hashCnt
          << " (literal activation)" << std::endl;
      for (; d_assertedHashes < static_cast<std::size_t>(hashCnt);
           ++d_assertedHashes)
      {
        hashPool[d_assertedHashes].assertToSolver(cs, d_useNativeXor);
      }
      auto& tm = ttc::getTermBuilder(cs);
      std::vector<cvc5::Term> assumptions;
      assumptions.reserve(static_cast<std::size_t>(hashCnt));
      for (std::int64_t i = 0; i < hashCnt; ++i)
      {
        const cvc5::Term& act = hashPool[i].activation();
        if (!act.isNull())
        {
          assumptions.push_back(tm.mkTerm(cvc5::Kind::NOT, {act}));
        }
      }
      attempt = d_counter.countWithAssumptions(assumptions, threshold);
    }
    else
    {
      std::vector<HashConstraint> active = activeHashes(hashCnt);
      Trace("pact") << "[pact] Evaluating hash count " << hashCnt << std::endl;
      // Native XOR clauses are added to CaDiCaL's Gaussian engine, which has no
      // per-scope retraction (an activation-guarded XOR cannot be deactivated --
      // forcing the activation only flips the parity). Reusing the solver across
      // galloping levels would therefore accumulate every level's hashes and
      // over-constrain the formula. Rebuild a fresh counting solver per level so
      // each count sees only its own hashes.
      if ((d_useNativeXor || d_useModpGj)
          && !std::getenv("TTC_NO_PERLEVEL_REBUILD"))
      {
        rebuildCountSolver();
      }

      // ApproxMC's reuse_models: pre-seed the count with the stored solutions
      // that still satisfy these `hashCnt` parities, so the SMT solver only has
      // to discover the ones we have not seen yet. A model saved with at least
      // `hashCnt` hashes active automatically satisfies this prefix; one saved
      // with fewer is checked against the parities (cheap term evaluation, no
      // SMT call). The fresh solver after a rebuild has no banning clauses, so
      // priming here is always sound. Reuse is skipped at level 0 (no hashes to
      // satisfy, and the base count already saturated).
      std::size_t primedCount = 0;
      if (d_reuseModels && hashCnt > 0)
      {
        std::vector<std::vector<cvc5::Term>> reusable;
        reusable.reserve(d_savedModels.size());
        for (const SavedModel& sm : d_savedModels)
        {
          ++curReuseChecked;
          if (sm.hashCount >= static_cast<std::size_t>(hashCnt)
              || modelSatisfiesHashes(sm.values, active))
          {
            reusable.push_back(sm.values);
            if (reusable.size() >= threshold)
            {
              break;
            }
          }
        }
        curReuseSat = reusable.size();
        primedCount = curReuseSat;
        d_counter.primeCache(active, std::move(reusable));
      }

      attempt = d_counter.count(active, threshold);

      // Record the newly enumerated solutions (those past the primed prefix) so
      // a later, lower hash count can reuse them in turn.
      if (d_reuseModels && hashCnt > 0)
      {
        const std::vector<std::vector<cvc5::Term>>& found =
            d_counter.lastModels();
        for (std::size_t i = primedCount; i < found.size(); ++i)
        {
          d_savedModels.push_back(
              {found[i], static_cast<std::size_t>(hashCnt)});
        }
      }
    }
    Trace("pact") << "[pact]   level " << hashCnt << " -> "
        << (attempt.has_value() ? std::to_string(*attempt)
                                : std::string(">=threshold"))
        << std::endl;
    const bool below = attempt.has_value();
    const std::int64_t numSols = below
                                     ? static_cast<std::int64_t>(*attempt)
                                     : static_cast<std::int64_t>(threshold) + 1;

    if (below)
    {
      // Found an exact count here. If one fewer hash saturated, the boundary is
      // right here and this is the answer.
      numExplored = lowerFib + totalMaxHashes - hashCnt;
      auto prev = thresholdSols.find(hashCnt - 1);
      if (hashCnt == 0 || (prev != thresholdSols.end() && prev->second))
      {
        prevMeasure = hashCnt;
        report(static_cast<std::size_t>(hashCnt), attempt, hashCnt, curReuseSat, curReuseChecked);
        Trace("pact") << "[pact] Winner at " << hashCnt << " with " << numSols
            << " models" << std::endl;
        return {static_cast<std::size_t>(hashCnt),
                static_cast<std::size_t>(numSols)};
      }

      thresholdSols[hashCnt] = false;
      solsForHash[hashCnt] = numSols;

      std::int64_t nextHash;
      if (iter > 0 && std::llabs(hashCnt - prevMeasure) <= 2)
      {
        // Close to last round's answer: step linearly (re-count) downward.
        upperFib = hashCnt;
        nextHash = hashCnt - 1;
      }
      else
      {
        if (hashPrev > hashCnt)
        {
          hashPrev = 0;
        }
        upperFib = hashCnt;
        if (hashPrev > lowerFib)
        {
          lowerFib = hashPrev;
        }
        nextHash = (upperFib + lowerFib) / 2;
      }
      report(static_cast<std::size_t>(hashCnt), attempt, nextHash, curReuseSat, curReuseChecked);
      hashPrev = curHashCnt;
      hashCnt = nextHash;
    }
    else
    {
      // Saturated the pivot here. If one more hash fell below, the boundary is
      // there and that count is the answer.
      numExplored = hashCnt + totalMaxHashes - upperFib;
      auto above = thresholdSols.find(hashCnt + 1);
      if (above != thresholdSols.end() && !above->second)
      {
        prevMeasure = hashCnt + 1;
        report(static_cast<std::size_t>(hashCnt), std::nullopt, hashCnt + 1, curReuseSat, curReuseChecked);
        Trace("pact") << "[pact] Winner at " << (hashCnt + 1) << " with "
            << solsForHash[hashCnt + 1] << " models" << std::endl;
        return {static_cast<std::size_t>(hashCnt + 1),
                static_cast<std::size_t>(solsForHash[hashCnt + 1])};
      }

      thresholdSols[hashCnt] = true;
      solsForHash[hashCnt] = static_cast<std::int64_t>(threshold) + 1;

      std::int64_t nextHash;
      if (iter > 0 && std::llabs(hashCnt - prevMeasure) < 2)
      {
        // Close to last round's answer: step linearly (re-count) upward.
        lowerFib = hashCnt;
        nextHash = hashCnt + 1;
      }
      else if (lowerFib + (hashCnt - lowerFib) * 2 >= upperFib - 1)
      {
        // Bracket is tight enough: switch to binary search.
        lowerFib = hashCnt;
        nextHash = (lowerFib + upperFib) / 2;
      }
      else
      {
        // Exponential (galloping) growth away from lowerFib.
        nextHash = lowerFib + (hashCnt - lowerFib) * 2;
        if (nextHash == hashCnt)
        {
          ++nextHash;
        }
      }
      report(static_cast<std::size_t>(hashCnt), std::nullopt, nextHash, curReuseSat, curReuseChecked);
      hashPrev = curHashCnt;
      hashCnt = nextHash;
    }
  }

  // Bracket fully explored without a sharp boundary (e.g. tiny variable sets):
  // fall back to the smallest hash count known to fall below the pivot.
  prevMeasure = upperFib;
  auto it = solsForHash.find(upperFib);
  std::int64_t cell = (it != solsForHash.end()) ? it->second : 1;
  return {static_cast<std::size_t>(upperFib), static_cast<std::size_t>(cell)};
}

Pact::MeasurementResult Pact::oneMeasurementAppmc7(
    std::size_t iter,
    std::int64_t& prevMeasure,
    std::vector<HashConstraint>& hashPool,
    const ReportFn& report)
{
  // ApproxMC7 searches for adjacent levels h and h+1 where h is non-empty and
  // h+1 is empty. A bounded count with threshold 1 returns nullopt on SAT
  // (non-empty) and exact 0 on UNSAT (empty).
  constexpr std::size_t existenceThreshold = 1;
  const std::int64_t totalMaxHashes = static_cast<std::int64_t>(d_maxHashes);
  std::int64_t lowerFib = 0;
  std::int64_t upperFib = totalMaxHashes + 1;
  std::int64_t hashCnt = prevMeasure;
  std::int64_t hashPrev = hashCnt;

  std::unordered_map<std::int64_t, bool> thresholdSols;
  std::unordered_map<std::int64_t, std::int64_t> solsForHash;
  thresholdSols[lowerFib] = true;
  solsForHash[lowerFib] = 1;
  thresholdSols[upperFib] = false;
  solsForHash[upperFib] = 0;

  auto activeHashes = [&](std::int64_t level) {
    while (static_cast<std::int64_t>(hashPool.size()) < level)
    {
      hashPool.push_back(generateHashConstraint());
    }
    return std::vector<HashConstraint>(hashPool.begin(),
                                       hashPool.begin() + level);
  };

  while (true)
  {
    if (hashCnt < 0)
    {
      hashCnt = 0;
    }
    if (hashCnt > totalMaxHashes)
    {
      hashCnt = totalMaxHashes;
    }

    const std::int64_t curHashCnt = hashCnt;
    std::optional<std::size_t> attempt;
    std::size_t curReuseSat = 0;
    std::size_t curReuseChecked = 0;

    if (d_xorActivationLiteral && d_useNativeXor)
    {
      activeHashes(hashCnt);
      cvc5::Solver& cs = d_countSolver ? *d_countSolver : d_solver;
      Trace("pact") << "[pact] Evaluating hash count " << hashCnt
          << " (ApproxMC7, literal activation)" << std::endl;
      for (; d_assertedHashes < static_cast<std::size_t>(hashCnt);
           ++d_assertedHashes)
      {
        hashPool[d_assertedHashes].assertToSolver(cs, d_useNativeXor);
      }
      auto& tm = ttc::getTermBuilder(cs);
      std::vector<cvc5::Term> assumptions;
      assumptions.reserve(static_cast<std::size_t>(hashCnt));
      for (std::int64_t i = 0; i < hashCnt; ++i)
      {
        const cvc5::Term& act = hashPool[i].activation();
        if (!act.isNull())
        {
          assumptions.push_back(tm.mkTerm(cvc5::Kind::NOT, {act}));
        }
      }
      attempt =
          d_counter.countWithAssumptions(assumptions, existenceThreshold);
    }
    else
    {
      std::vector<HashConstraint> active = activeHashes(hashCnt);
      Trace("pact") << "[pact] Evaluating hash count " << hashCnt
          << " (ApproxMC7)" << std::endl;
      if ((d_useNativeXor || d_useModpGj)
          && !std::getenv("TTC_NO_PERLEVEL_REBUILD"))
      {
        rebuildCountSolver();
      }

      std::size_t primedCount = 0;
      if (d_reuseModels && hashCnt > 0)
      {
        std::vector<std::vector<cvc5::Term>> reusable;
        reusable.reserve(d_savedModels.size());
        for (const SavedModel& sm : d_savedModels)
        {
          ++curReuseChecked;
          if (sm.hashCount >= static_cast<std::size_t>(hashCnt)
              || modelSatisfiesHashes(sm.values, active))
          {
            reusable.push_back(sm.values);
            break;
          }
        }
        curReuseSat = reusable.size();
        primedCount = curReuseSat;
        d_counter.primeCache(active, std::move(reusable));
      }

      attempt = d_counter.count(active, existenceThreshold);

      if (d_reuseModels && hashCnt > 0)
      {
        const std::vector<std::vector<cvc5::Term>>& found =
            d_counter.lastModels();
        for (std::size_t i = primedCount; i < found.size(); ++i)
        {
          d_savedModels.push_back(
              {found[i], static_cast<std::size_t>(hashCnt)});
        }
      }
    }

    const bool empty = attempt.has_value() && *attempt == 0;
    Trace("pact") << "[pact]   level " << hashCnt << " -> "
        << (empty ? std::string("0") : std::string("non-empty"))
        << std::endl;

    if (empty)
    {
      if (hashCnt == 0)
      {
        prevMeasure = 0;
        report(0, std::size_t{0}, 0, curReuseSat, curReuseChecked);
        Trace("pact") << "[pact] ApproxMC7 found UNSAT at level 0"
            << std::endl;
        return {0, 0};
      }
      auto prev = thresholdSols.find(hashCnt - 1);
      if (prev != thresholdSols.end() && prev->second)
      {
        prevMeasure = hashCnt - 1;
        report(static_cast<std::size_t>(hashCnt),
               std::size_t{0},
               hashCnt - 1,
               curReuseSat,
               curReuseChecked);
        Trace("pact") << "[pact] ApproxMC7 winner at " << (hashCnt - 1)
            << std::endl;
        return {static_cast<std::size_t>(hashCnt - 1), 1};
      }

      thresholdSols[hashCnt] = false;
      solsForHash[hashCnt] = 0;

      std::int64_t nextHash;
      if (iter > 0 && std::llabs(hashCnt - prevMeasure) <= 2)
      {
        upperFib = hashCnt;
        nextHash = hashCnt - 1;
      }
      else
      {
        if (hashPrev > hashCnt)
        {
          hashPrev = 0;
        }
        upperFib = hashCnt;
        if (hashPrev > lowerFib)
        {
          lowerFib = hashPrev;
        }
        nextHash = (upperFib + lowerFib) / 2;
      }
      report(static_cast<std::size_t>(hashCnt),
             std::size_t{0},
             nextHash,
             curReuseSat,
             curReuseChecked);
      hashPrev = curHashCnt;
      hashCnt = nextHash;
    }
    else
    {
      auto above = thresholdSols.find(hashCnt + 1);
      if (above != thresholdSols.end() && !above->second)
      {
        prevMeasure = hashCnt;
        report(static_cast<std::size_t>(hashCnt),
               std::nullopt,
               hashCnt,
               curReuseSat,
               curReuseChecked);
        Trace("pact") << "[pact] ApproxMC7 winner at " << hashCnt
            << std::endl;
        return {static_cast<std::size_t>(hashCnt), 1};
      }

      thresholdSols[hashCnt] = true;
      solsForHash[hashCnt] = 1;

      std::int64_t nextHash;
      if (iter > 0 && std::llabs(hashCnt - prevMeasure) < 2)
      {
        lowerFib = hashCnt;
        nextHash = hashCnt + 1;
      }
      else if (lowerFib + (hashCnt - lowerFib) * 2 >= upperFib - 1)
      {
        lowerFib = hashCnt;
        nextHash = (lowerFib + upperFib) / 2;
      }
      else
      {
        nextHash = lowerFib + (hashCnt - lowerFib) * 2;
        if (nextHash == hashCnt)
        {
          ++nextHash;
        }
      }
      report(static_cast<std::size_t>(hashCnt),
             std::nullopt,
             nextHash,
             curReuseSat,
             curReuseChecked);
      hashPrev = curHashCnt;
      hashCnt = nextHash;
    }
  }
}

bool Pact::modelSatisfiesHashes(
    const std::vector<cvc5::Term>& model,
    const std::vector<HashConstraint>& hashes)
{
  // Evaluate each parity under the saved projection assignment without touching
  // the SAT solver. Every term in a hash's XOR clauses references only
  // projection variables (a Boolean projection var, or a (= (extract b b) #b1)
  // bit literal of a projection bit-vector), so substituting the model values
  // and simplifying folds it to a Boolean constant. The hash holds iff the XOR
  // of those constants equals the clause's right-hand side.
  cvc5::Solver& s = d_countSolver ? *d_countSolver : d_solver;
  for (const HashConstraint& hash : hashes)
  {
    // Word-level (Prime/Lemire) hashes carry no parity literals -- their clause
    // has empty `terms` and the whole constraint lives in the fallback equality.
    // Evaluate that term under the saved assignment instead of XORing literals.
    bool arithmetic = false;
    for (const XorClause& clause : hash.xorClauses())
    {
      // Word-level Prime/Lemire have empty parity terms; PrimeGj has parity-like
      // bit literals but is really a mod-p row (modulus != 0). Both are checked
      // by evaluating the (arithmetic) fallback term, not by XORing literals.
      if (clause.terms.empty() || clause.modulus != 0)
      {
        arithmetic = true;
        break;
      }
    }
    if (arithmetic)
    {
      cvc5::Term value =
          s.simplify(hash.fallback().substitute(d_projectionVars, model));
      if (!(value.isBooleanValue() && value.getBooleanValue()))
      {
        return false;
      }
      continue;
    }
    for (const XorClause& clause : hash.xorClauses())
    {
      bool parity = false;
      for (const cvc5::Term& term : clause.terms)
      {
        cvc5::Term value =
            s.simplify(term.substitute(d_projectionVars, model));
        if (value.isBooleanValue() && value.getBooleanValue())
        {
          parity = !parity;
        }
      }
      if (parity != clause.rhs)
      {
        return false;
      }
    }
  }
  return true;
}
