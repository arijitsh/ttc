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
           bool bvPact)
    : d_solver(solver),
      d_rng(static_cast<std::mt19937::result_type>(seed)),
      d_counter(solver),
      d_boolHash(solver),
      d_bvHash(solver)
{
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
      "sat-solver",
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
}

Pact::Parameters Pact::getParameters() const
{
  const double epsilon = 0.8;
  const double delta = 2.5;
  std::size_t pivot = static_cast<std::size_t>(
      std::ceil(2.0 * std::ceil(4.49 * std::pow(1.0 + 1.0 / epsilon, 2.0))));
  if (pivot == 0)
  {
    pivot = 1;
  }
  if (const char* e = std::getenv("TTC_PIVOT"))
  {
    pivot = static_cast<std::size_t>(std::atoll(e));
  }
  std::size_t iterations = static_cast<std::size_t>(
      std::ceil(25.0 * std::log(3.0 / delta)));
  if (iterations == 0)
  {
    iterations = 1;
  }
  return {pivot, iterations};
}

std::uint64_t Pact::median(std::vector<std::uint64_t>& values) const
{
  if (values.empty())
  {
    return 0;
  }
  std::size_t mid = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + mid, values.end());
  return values[mid];
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

  return HashConstraint(fallback, std::move(clauses));
}

std::uint64_t Pact::count()
{
  Parameters params = getParameters();
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
                    std::int64_t nextHash) {
    if (hashCount > maxHashesUsed)
    {
      maxHashesUsed = hashCount;
    }
    if (Log.getVerbosity() == 0)
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
      countStream << ">= " << params.threshold;
    }
    std::ostringstream line;
    line << std::fixed << std::setprecision(3);
    line << "c " << std::setw(8) << elapsed;
    line << "  " << std::setw(6) << currentRound;
    line << "  " << std::setw(6) << hashCount;
    line << "  " << std::setw(10) << countStream.str();
    line << "  " << std::setw(9) << nextHash;
    std::cout << line.str() << std::endl;
  };

  // Base count with no hashing. If the formula has fewer than `threshold`
  // models in total the count is exact and we are done -- this mirrors ApproxMC
  // "counting without XORs".
  rebuildCountSolver();
  std::vector<HashConstraint> empty;
  std::optional<std::size_t> base = d_counter.count(empty, params.threshold);
  report(0, base, 0);
  if (base.has_value())
  {
    Trace("pact") << "[pact] Base model count succeeded without hashing: "
        << *base << std::endl;
    std::cout << "c max hashes: " << maxHashesUsed << std::endl;
    return static_cast<std::uint64_t>(*base);
  }
  Trace("pact")
      << "[pact] Base count reached saturation; introducing random hashes"
      << std::endl;

  std::vector<std::uint64_t> estimates;
  estimates.reserve(params.iterations);

  // Carries the winning hash count between measurements (ApproxMC's
  // prev_measure) so later rounds start their galloping search near the answer.
  std::int64_t prevMeasure = 0;

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
    // deterministic.
    std::vector<HashConstraint> hashPool;
    MeasurementResult m =
        oneMeasurement(iter, prevMeasure, hashPool, params.threshold, report);

    std::uint64_t approx = static_cast<std::uint64_t>(m.cellCount);
    approx <<= m.hashCount;
    Trace("pact") << "[pact] Iteration " << (iter + 1) << " estimate: "
        << approx << " (" << m.cellCount << " << " << m.hashCount << ")"
        << std::endl;
    estimates.push_back(approx);

    if (Log.getVerbosity() != 0)
    {
      double elapsed = Log.elapsed() - startTime;
      std::ostringstream line;
      line << std::fixed << std::setprecision(3);
      line << "c " << std::setw(8) << elapsed;
      line << "  " << std::setw(6) << currentRound;
      line << "  " << std::setw(6) << m.hashCount;
      line << "  " << std::setw(10) << m.cellCount;
      line << "  " << std::setw(9) << "winner";
      line << "  " << std::setw(16) << approx;
      std::cout << line.str() << std::endl;
      lastProgress = elapsed;
    }
  }

  std::uint64_t result = median(estimates);
  Trace("pact") << "[pact] Median estimate across iterations: " << result
      << std::endl;
  std::cout << "c max hashes: " << maxHashesUsed << std::endl;
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
    std::vector<HashConstraint> active = activeHashes(hashCnt);
    Trace("pact") << "[pact] Evaluating hash count " << hashCnt << std::endl;
    // Native XOR clauses are added to CaDiCaL's Gaussian engine, which has no
    // per-scope retraction (an activation-guarded XOR cannot be deactivated --
    // forcing the activation only flips the parity). Reusing the solver across
    // galloping levels would therefore accumulate every level's hashes and
    // over-constrain the formula. Rebuild a fresh counting solver per level so
    // each count sees only its own hashes.
    if (d_useNativeXor && !std::getenv("TTC_NO_PERLEVEL_REBUILD"))
    {
      rebuildCountSolver();
    }
    std::optional<std::size_t> attempt = d_counter.count(active, threshold);
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
        report(static_cast<std::size_t>(hashCnt), attempt, hashCnt);
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
      report(static_cast<std::size_t>(hashCnt), attempt, nextHash);
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
        report(static_cast<std::size_t>(hashCnt), std::nullopt, hashCnt + 1);
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
      report(static_cast<std::size_t>(hashCnt), std::nullopt, nextHash);
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

