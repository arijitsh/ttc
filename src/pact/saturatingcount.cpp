#include "pact/saturatingcount.hpp"

#include "features.hpp"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <utility>
#include <set>
#include <vector>

#include "logger.hpp"

SaturatingCounter::SaturatingCounter(
    cvc5::Solver& solver, const std::vector<cvc5::Term>* projectionVars)
    : d_solver(&solver), d_projectionVars(projectionVars)
{}

void SaturatingCounter::setProjectionVars(
    const std::vector<cvc5::Term>* projectionVars)
{
  d_projectionVars = projectionVars;
}

void SaturatingCounter::resetStatistics()
{
  d_smtCalls = 0;
}

void SaturatingCounter::resetCache()
{
  d_cachedConstraints.clear();
  d_cachedModels.clear();
}

std::optional<std::size_t> SaturatingCounter::count(
    const std::vector<HashConstraint>& additionalConstraints,
    std::size_t threshold)
{
  if (!d_projectionVars)
  {
    throw std::logic_error("Projection variables not initialised for counter");
  }
  Trace("saturating") << "Counting with " << additionalConstraints.size()
      << " additional constraints and threshold " << threshold << std::endl;
  if (threshold == 0)
  {
    Trace("saturating") << "Threshold is zero; returning 0 immediately"
        << std::endl;
    return std::size_t{0};
  }

  if (d_projectionVars->empty())
  {
    Trace("saturating") << "No projection variables; single satisfiability"
        << " check" << std::endl;
    d_solver->push();
    for (const HashConstraint& constraint : additionalConstraints)
    {
      constraint.assertToSolver(*d_solver, d_useNativeXor);
    }
    ++d_smtCalls;
    cvc5::Result res = d_solver->checkSat();
    d_solver->pop();
    if (res.isSat())
    {
      Trace("saturating") << "Formula is SAT without projections" << std::endl;
      return std::size_t{1};
    }
    if (res.isUnsat())
    {
      Trace("saturating") << "Formula is UNSAT without projections"
          << std::endl;
      return std::size_t{0};
    }
    Trace("saturating") << "Solver returned unknown without projections"
        << std::endl;
    return std::nullopt;
  }

  Trace("saturating") << "Enumerating assignments for "
      << d_projectionVars->size() << " projection variables" << std::endl;
  const bool projectionsAreBoolean = std::all_of(
      d_projectionVars->cbegin(),
      d_projectionVars->cend(),
      [](const cvc5::Term& term) { return term.getSort().isBoolean(); });
  auto& tm = ttc::getTermBuilder(*d_solver);
  d_solver->push();
  for (const HashConstraint& constraint : additionalConstraints)
  {
    constraint.assertToSolver(*d_solver, d_useNativeXor);
  }

  std::optional<std::size_t> result;
  std::size_t modelCount = 0;
  std::vector<std::vector<cvc5::Term>> currentModels;
  currentModels.reserve(threshold);
  std::set<std::string> dupSeen;  // within-count duplicate detection (debug)

  auto buildBlocking = [&](const std::vector<cvc5::Term>& modelValues,
                           std::size_t index) {
    bool useBooleanBlocking = projectionsAreBoolean;
    if (useBooleanBlocking)
    {
      for (const cvc5::Term& value : modelValues)
      {
        if (!value.isBooleanValue())
        {
          useBooleanBlocking = false;
          break;
        }
      }
    }
    if (useBooleanBlocking)
    {
      std::vector<cvc5::Term> clause;
      clause.reserve(modelValues.size());
      for (std::size_t i = 0; i < modelValues.size(); ++i)
      {
        const cvc5::Term& var = (*d_projectionVars)[i];
        const cvc5::Term& value = modelValues[i];
        clause.push_back(value.getBooleanValue()
                             ? tm.mkTerm(cvc5::Kind::NOT, {var})
                             : var);
      }
      Trace("saturating") << "Blocking model " << index << " with "
          << clause.size() << " boolean literals" << std::endl;
      if (clause.size() == 1)
      {
        return clause.front();
      }
      return tm.mkTerm(cvc5::Kind::OR, clause);
    }
    std::vector<cvc5::Term> equalities;
    equalities.reserve(modelValues.size());
    for (std::size_t i = 0; i < modelValues.size(); ++i)
    {
      equalities.push_back(tm.mkTerm(
          cvc5::Kind::EQUAL, {(*d_projectionVars)[i], modelValues[i]}));
    }
    Trace("saturating") << "Blocking model " << index << " with "
        << equalities.size() << " equalities" << std::endl;
    cvc5::Term blocking;
    if (equalities.size() == 1)
    {
      blocking = equalities.front();
    }
    else
    {
      blocking = tm.mkTerm(cvc5::Kind::AND, equalities);
    }
    return tm.mkTerm(cvc5::Kind::NOT, {blocking});
  };

  std::size_t prefix = 0;
  while (prefix < d_cachedConstraints.size()
         && prefix < additionalConstraints.size()
         && d_cachedConstraints[prefix] == additionalConstraints[prefix])
  {
    ++prefix;
  }

  bool reuseCached = false;
  bool filterCached = false;
  if (!d_cachedConstraints.empty() || !d_cachedModels.empty())
  {
    if (prefix == d_cachedConstraints.size()
        && prefix == additionalConstraints.size())
    {
      reuseCached = true;
    }
    else if (prefix == d_cachedConstraints.size()
             && additionalConstraints.size() > d_cachedConstraints.size())
    {
      reuseCached = true;
      filterCached = true;
    }
    else if (prefix == additionalConstraints.size()
             && d_cachedConstraints.size() > additionalConstraints.size())
    {
      reuseCached = true;
    }
    else
    {
      d_cachedConstraints.clear();
      d_cachedModels.clear();
    }
  }

  if (reuseCached)
  {
    if (filterCached)
    {
      std::vector<std::vector<cvc5::Term>> survivors;
      survivors.reserve(d_cachedModels.size());
      for (const auto& model : d_cachedModels)
      {
        d_solver->push();
        for (std::size_t i = 0; i < model.size(); ++i)
        {
          cvc5::Term eq = tm.mkTerm(
              cvc5::Kind::EQUAL, {(*d_projectionVars)[i], model[i]});
          d_solver->assertFormula(eq);
        }
        ++d_smtCalls;
        cvc5::Result res = d_solver->checkSat();
        d_solver->pop();
        if (!res.isSat())
        {
          continue;
        }
        ++modelCount;
        survivors.push_back(model);
        currentModels.push_back(model);
        cvc5::Term blocking = buildBlocking(model, modelCount);
        d_solver->assertFormula(blocking);
        if (modelCount >= threshold)
        {
          Trace("saturating") << "Cached models reached threshold "
              << threshold << std::endl;
          result = std::nullopt;
          break;
        }
      }
      d_cachedModels = std::move(survivors);
    }
    else
    {
      for (const auto& model : d_cachedModels)
      {
        ++modelCount;
        currentModels.push_back(model);
        cvc5::Term blocking = buildBlocking(model, modelCount);
        d_solver->assertFormula(blocking);
        if (modelCount >= threshold)
        {
          Trace("saturating") << "Cached models reached threshold "
              << threshold << std::endl;
          result = std::nullopt;
          break;
        }
      }
    }
  }

  if (!result.has_value() && modelCount < threshold)
  {
    while (true)
    {
      std::size_t callIndex = modelCount + 1;
      Trace("saturating") << "checkSat call " << callIndex << std::endl;
      ++d_smtCalls;
      cvc5::Result res = d_solver->checkSat();
      if (res.isUnsat())
      {
        Trace("saturating") << "UNSAT after " << modelCount
            << " models" << std::endl;
        result = modelCount;
        break;
      }
      if (!res.isSat())
      {
        Trace("saturating") << "Solver returned unknown after "
            << modelCount << " models" << std::endl;
        result = std::nullopt;
        break;
      }

      ++modelCount;
      Trace("saturating") << "Recorded model " << modelCount << std::endl;
      if (modelCount >= threshold)
      {
        Trace("saturating") << "Reached threshold " << threshold
            << "; stopping enumeration" << std::endl;
        result = std::nullopt;
        break;
      }

      std::vector<cvc5::Term> modelValues;
      modelValues.reserve(d_projectionVars->size());
      for (const cvc5::Term& var : *d_projectionVars)
      {
        modelValues.push_back(d_solver->getValue(var));
      }
      if (std::getenv("TTC_DUP_DEBUG"))
      {
        std::string key;
        for (const auto& v : modelValues) key += v.toString() + ",";
        if (!dupSeen.insert(key).second)
          std::cerr << "[dup] WITHIN-COUNT duplicate at model " << modelCount
                    << " (distinct=" << dupSeen.size() << ")" << std::endl;
      }
      if (std::getenv("TTC_XOR_VERIFY"))
      {
        for (const HashConstraint& hc : additionalConstraints)
        {
          for (const XorClause& xc : hc.xorClauses())
          {
            if (xc.terms.empty()) continue;
            bool parity = false;
            for (const cvc5::Term& t : xc.terms)
            {
              cvc5::Term val = d_solver->getValue(t);
              if (val.isBooleanValue() && val.getBooleanValue())
                parity = !parity;
            }
            if (parity != xc.rhs)
            {
              std::cerr << "[xorverify] model " << modelCount
                        << " VIOLATES a hash XOR (parity=" << parity
                        << " rhs=" << xc.rhs << ")" << std::endl;
              break;
            }
          }
        }
      }
      currentModels.push_back(modelValues);
      cvc5::Term blockingConstraint = buildBlocking(modelValues, modelCount);
      d_solver->assertFormula(blockingConstraint);
    }
  }

  d_solver->pop();
  if (std::getenv("TTC_COUNT_DEBUG"))
    std::cerr << "[count] hashes=" << additionalConstraints.size()
              << " modelCount=" << modelCount << " result="
              << (result.has_value() ? std::to_string(*result)
                                      : std::string("saturated"))
              << std::endl;
  d_cachedConstraints.assign(additionalConstraints.begin(),
                             additionalConstraints.end());
  d_cachedModels = currentModels;
  if (result.has_value())
  {
    Trace("saturating") << "Enumerated " << *result
        << " assignments below the threshold" << std::endl;
  }
  else
  {
    Trace("saturating") << "Enumeration terminated without an exact count"
        << std::endl;
  }
  return result;
}

cvc5::Term SaturatingCounter::buildBlockingClause(
    const std::vector<cvc5::Term>& modelValues, bool projectionsAreBoolean) const
{
  auto& tm = ttc::getTermBuilder(*d_solver);
  bool useBooleanBlocking = projectionsAreBoolean;
  if (useBooleanBlocking)
  {
    for (const cvc5::Term& value : modelValues)
    {
      if (!value.isBooleanValue())
      {
        useBooleanBlocking = false;
        break;
      }
    }
  }
  if (useBooleanBlocking)
  {
    std::vector<cvc5::Term> clause;
    clause.reserve(modelValues.size());
    for (std::size_t i = 0; i < modelValues.size(); ++i)
    {
      const cvc5::Term& var = (*d_projectionVars)[i];
      clause.push_back(modelValues[i].getBooleanValue()
                           ? tm.mkTerm(cvc5::Kind::NOT, {var})
                           : var);
    }
    if (clause.size() == 1)
    {
      return clause.front();
    }
    return tm.mkTerm(cvc5::Kind::OR, clause);
  }
  std::vector<cvc5::Term> equalities;
  equalities.reserve(modelValues.size());
  for (std::size_t i = 0; i < modelValues.size(); ++i)
  {
    equalities.push_back(
        tm.mkTerm(cvc5::Kind::EQUAL, {(*d_projectionVars)[i], modelValues[i]}));
  }
  cvc5::Term blocking = equalities.size() == 1
                            ? equalities.front()
                            : tm.mkTerm(cvc5::Kind::AND, equalities);
  return tm.mkTerm(cvc5::Kind::NOT, {blocking});
}

std::optional<std::size_t> SaturatingCounter::countWithAssumptions(
    const std::vector<cvc5::Term>& assumptions, std::size_t threshold)
{
  if (!d_projectionVars)
  {
    throw std::logic_error("Projection variables not initialised for counter");
  }
  if (threshold == 0)
  {
    return std::size_t{0};
  }

  // The hash parities are already on the solver; the active ones are selected by
  // 'assumptions' (each an asserted-false indicator). A fresh measurement leaves
  // no cache to reuse here, so this path keeps no model cache.
  auto sat = [&]() {
    ++d_smtCalls;
    return assumptions.empty() ? d_solver->checkSat()
                               : d_solver->checkSatAssuming(assumptions);
  };

  if (d_projectionVars->empty())
  {
    cvc5::Result res = sat();
    if (res.isSat()) return std::size_t{1};
    if (res.isUnsat()) return std::size_t{0};
    return std::nullopt;
  }

  const bool projectionsAreBoolean = std::all_of(
      d_projectionVars->cbegin(), d_projectionVars->cend(),
      [](const cvc5::Term& term) { return term.getSort().isBoolean(); });

  d_solver->push();
  std::optional<std::size_t> result;
  std::size_t modelCount = 0;
  while (true)
  {
    cvc5::Result res = sat();
    if (res.isUnsat())
    {
      result = modelCount;
      break;
    }
    if (!res.isSat())
    {
      result = std::nullopt;
      break;
    }
    ++modelCount;
    if (modelCount >= threshold)
    {
      result = std::nullopt;
      break;
    }
    std::vector<cvc5::Term> modelValues;
    modelValues.reserve(d_projectionVars->size());
    for (const cvc5::Term& var : *d_projectionVars)
    {
      modelValues.push_back(d_solver->getValue(var));
    }
    d_solver->assertFormula(
        buildBlockingClause(modelValues, projectionsAreBoolean));
  }
  d_solver->pop();
  return result;
}

