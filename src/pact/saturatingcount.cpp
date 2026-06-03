#include "pact/saturatingcount.hpp"

#include "features.hpp"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "logger.hpp"

SaturatingCounter::SaturatingCounter(
    cvc5::Solver& solver, const std::vector<cvc5::Term>* projectionVars)
    : d_solver(solver), d_projectionVars(projectionVars)
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
    d_solver.push();
    for (const HashConstraint& constraint : additionalConstraints)
    {
      constraint.assertToSolver(d_solver, d_useNativeXor);
    }
    ++d_smtCalls;
    cvc5::Result res = d_solver.checkSat();
    d_solver.pop();
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
  auto& tm = ttc::getTermBuilder(d_solver);
  if (std::getenv("TTC_NO_CACHE") != nullptr)
  {
    d_cachedConstraints.clear();
    d_cachedModels.clear();
  }
  d_solver.push();
  for (const HashConstraint& constraint : additionalConstraints)
  {
    constraint.assertToSolver(d_solver, d_useNativeXor);
  }

  std::optional<std::size_t> result;
  std::size_t modelCount = 0;
  std::vector<std::vector<cvc5::Term>> currentModels;
  currentModels.reserve(threshold);

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
        d_solver.push();
        for (std::size_t i = 0; i < model.size(); ++i)
        {
          cvc5::Term eq = tm.mkTerm(
              cvc5::Kind::EQUAL, {(*d_projectionVars)[i], model[i]});
          d_solver.assertFormula(eq);
        }
        ++d_smtCalls;
        cvc5::Result res = d_solver.checkSat();
        d_solver.pop();
        if (!res.isSat())
        {
          continue;
        }
        ++modelCount;
        survivors.push_back(model);
        currentModels.push_back(model);
        cvc5::Term blocking = buildBlocking(model, modelCount);
        d_solver.assertFormula(blocking);
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
        d_solver.assertFormula(blocking);
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
      cvc5::Result res = d_solver.checkSat();
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
        modelValues.push_back(d_solver.getValue(var));
      }
      currentModels.push_back(modelValues);
      cvc5::Term blockingConstraint = buildBlocking(modelValues, modelCount);
      d_solver.assertFormula(blockingConstraint);
    }
  }

  d_solver.pop();
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

