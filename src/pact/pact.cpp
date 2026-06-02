#include "pact/pact.hpp"

#include <algorithm>
#include <cmath>
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
           bool useNativeXor)
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

  d_counter.setProjectionVars(&d_projectionVars);
  d_counter.setUseNativeXor(useNativeXor);
  d_boolHash.setVariables(d_booleanVars, d_booleanBvVars);
  d_bvHash.setVariables(d_bitvectorVars);
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
  auto getNextHashDelta = [&](std::size_t currentHashes) -> std::size_t {
    if (!d_nextIndex.needsRefinement())
    {
      return 0;
    }
    std::optional<std::size_t> nextLevel = d_nextIndex.previewNextCandidate();
    if (!nextLevel.has_value())
    {
      return 0;
    }
    if (*nextLevel > currentHashes)
    {
      return *nextLevel - currentHashes;
    }
    return 0;
  };
  auto printProgress = [&](std::size_t round,
                           std::size_t hashCount,
                           const std::optional<std::size_t>& countResult,
                           std::size_t nextHashDelta,
                           bool force,
                           const std::optional<std::uint64_t>& roundEstimate) -> bool {
    if (Log.getVerbosity() == 0)
    {
      return false;
    }
    double elapsed = Log.elapsed() - startTime;
    bool shouldPrint = force || countResult.has_value();
    if (!shouldPrint && (elapsed - lastProgress < progress_interval))
    {
      return false;
    }
    lastProgress = elapsed;
    std::ostringstream countStream;
    if (countResult.has_value())
    {
      countStream << *countResult;
    }
    else
    {
      countStream << ">= " << params.threshold;
    }
    std::ostringstream line;
    line << std::fixed << std::setprecision(3);
    line << "c " << std::setw(8) << elapsed;
    line << "  " << std::setw(6) << round;
    line << "  " << std::setw(8) << hashCount;
    line << "  " << std::setw(16) << countStream.str();
    line << "  " << std::setw(9) << nextHashDelta;
    std::string estimateStr;
    if (roundEstimate.has_value())
    {
      std::ostringstream estStream;
      estStream << *roundEstimate;
      estimateStr = estStream.str();
    }
    line << "  " << std::setw(16) << estimateStr;
    std::cout << line.str() << std::endl;
    return true;
  };

  std::vector<HashConstraint> empty;
  std::size_t initialLevel = d_nextIndex.initial();
  std::optional<std::size_t> base = d_counter.count(empty, params.threshold);
  std::optional<std::uint64_t> baseEstimate;
  if (base.has_value())
  {
    baseEstimate = static_cast<std::uint64_t>(*base);
  }
  printProgress(0, 0, base, getNextHashDelta(0), true, baseEstimate);
  if (base.has_value())
  {
    Trace("pact") << "[pact] Base model count succeeded without hashing: " << *base
        << std::endl;
    return static_cast<std::uint64_t>(*base);
  }
  Trace("pact")
      << "[pact] Base count reached saturation; introducing random hashes"
      << std::endl;
  d_nextIndex.updateOnSaturation(initialLevel);

  std::vector<std::uint64_t> estimates;
  estimates.reserve(params.iterations);

  for (std::size_t iter = 0; iter < params.iterations; ++iter)
  {
    Trace("pact") << "[pact] Iteration " << (iter + 1) << " of "
        << params.iterations << std::endl;
    d_nextIndex.startRound();
    std::vector<HashConstraint> hashPool;
    std::vector<HashConstraint> activeHashes;
    std::unordered_map<std::size_t, std::optional<std::size_t>> levelResults;
    std::optional<std::size_t> bestCount;
    std::size_t bestLevel = 0;

    auto prepareHashes = [&](std::size_t level) {
      if (level > hashPool.size())
      {
        while (hashPool.size() < level)
        {
          HashConstraint h = generateHashConstraint();
          hashPool.push_back(h);
        }
      }
      activeHashes.assign(hashPool.begin(), hashPool.begin() + level);
    };

    while (true)
    {
      std::size_t level = d_nextIndex.nextCandidate();
      auto cacheIt = levelResults.find(level);
      if (cacheIt != levelResults.end())
      {
        Trace("pact") << "[pact] Reusing cached result for hash level " << level
            << std::endl;
        const std::optional<std::size_t>& cached = cacheIt->second;
        if (cached.has_value())
        {
          bestCount = cached;
          bestLevel = level;
          d_nextIndex.updateOnSatisfying(level);
          if (!d_nextIndex.needsRefinement())
          {
            break;
          }
        }
        else
        {
          d_nextIndex.updateOnSaturation(level);
        }
        continue;
      }

      prepareHashes(level);
      Trace("pact") << "[pact] Evaluating hash level " << level << std::endl;
      std::optional<std::size_t> attempt =
          d_counter.count(activeHashes, params.threshold);
      levelResults.emplace(level, attempt);
      std::optional<std::uint64_t> roundEstimate;
      if (attempt.has_value())
      {
        bestCount = attempt;
        bestLevel = level;
        Trace("pact") << "[pact] Hash level " << level << " yielded " << *attempt
            << " models" << std::endl;
        d_nextIndex.updateOnSatisfying(level);
        bool done = !d_nextIndex.needsRefinement();
        if (done)
        {
          roundEstimate = static_cast<std::uint64_t>(*attempt) << level;
        }
        bool forcePrint = roundEstimate.has_value();
        printProgress(iter + 1,
                      activeHashes.size(),
                      attempt,
                      getNextHashDelta(activeHashes.size()),
                      forcePrint,
                      roundEstimate);
        if (done)
        {
          break;
        }
      }
      else
      {
        Trace("pact") << "[pact] Hash level " << level
            << " saturated pivot " << params.threshold << std::endl;
        d_nextIndex.updateOnSaturation(level);
        bool done = bestCount.has_value() && !d_nextIndex.needsRefinement();
        if (done)
        {
          roundEstimate = static_cast<std::uint64_t>(*bestCount) << bestLevel;
        }
        bool forcePrint = true;
        printProgress(iter + 1,
                      activeHashes.size(),
                      attempt,
                      getNextHashDelta(activeHashes.size()),
                      forcePrint,
                      roundEstimate);
        if (done)
        {
          break;
        }
      }
    }

    std::uint64_t approx = static_cast<std::uint64_t>(*bestCount);
    approx <<= bestLevel;
    Trace("pact") << "[pact] Iteration " << (iter + 1)
        << " estimate: " << approx << std::endl;
    estimates.push_back(approx);
    d_nextIndex.finishRound();
    if (bestCount.has_value())
    {
      d_nextIndex.updateOnSatisfying(bestLevel);
    }
  }

  std::uint64_t result = median(estimates);
  Trace("pact") << "[pact] Median estimate across iterations: " << result
      << std::endl;
  return result;
}

