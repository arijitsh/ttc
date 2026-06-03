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

  // Emits one progress row per saturating-count evaluation, throttled by
  // progress_interval. `nextHash` is where the galloping search will jump next.
  auto report = [&](std::size_t hashCount,
                    const std::optional<std::size_t>& sols,
                    std::int64_t nextHash) {
    if (Log.getVerbosity() == 0)
    {
      return;
    }
    double elapsed = Log.elapsed() - startTime;
    if (elapsed - lastProgress < progress_interval)
    {
      return;
    }
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
    line << "  " << std::setw(8) << hashCount;
    line << "  " << std::setw(16) << countStream.str();
    line << "  " << std::setw(9) << nextHash;
    std::cout << line.str() << std::endl;
  };

  // Base count with no hashing. If the formula has fewer than `threshold`
  // models in total the count is exact and we are done -- this mirrors ApproxMC
  // "counting without XORs".
  std::vector<HashConstraint> empty;
  std::optional<std::size_t> base = d_counter.count(empty, params.threshold);
  report(0, base, 0);
  if (base.has_value())
  {
    Trace("pact") << "[pact] Base model count succeeded without hashing: "
        << *base << std::endl;
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
      line << "  " << std::setw(8) << m.hashCount;
      line << "  " << std::setw(16) << m.cellCount;
      line << "  " << std::setw(9) << "winner";
      line << "  " << std::setw(16) << approx;
      std::cout << line.str() << std::endl;
      lastProgress = elapsed;
    }
  }

  std::uint64_t result = median(estimates);
  Trace("pact") << "[pact] Median estimate across iterations: " << result
      << std::endl;
  if (std::getenv("TTC_PACT_STATS") != nullptr)
  {
    try
    {
      std::string stats = d_solver.getStatistics().toString();
      std::istringstream iss(stats);
      std::string line;
      while (std::getline(iss, line))
      {
        if (line.find("cryptominisat") != std::string::npos
            || line.find("BVSolverBitblast") != std::string::npos
            || line.find("cadical") != std::string::npos
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
    std::optional<std::size_t> attempt = d_counter.count(active, threshold);
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

