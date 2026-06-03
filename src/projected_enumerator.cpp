#include "projected_enumerator.hpp"

#include <stdexcept>

#include "features.hpp"

ProjectedEnumerator::ProjectedEnumerator(
    cvc5::Solver& solver, const std::vector<cvc5::Term>& projectionVars)
    : d_solver(solver), d_projectionVars(projectionVars), d_smtCalls(0)
{
}

std::uint64_t ProjectedEnumerator::count()
{
  d_smtCalls = 0;
  d_solver.push();
  std::uint64_t total = 0;
  try
  {
    if (d_projectionVars.empty())
    {
      cvc5::Result res = d_solver.checkSat();
      ++d_smtCalls;
      if (res.isSat())
      {
        total = 1;
      }
      else if (res.isUnsat())
      {
        total = 0;
      }
      else
      {
        throw std::runtime_error("solver returned unknown during enumeration");
      }
      d_solver.pop();
      return total;
    }

    auto& tm = ttc::getTermBuilder(d_solver);
    while (true)
    {
      cvc5::Result res = d_solver.checkSat();
      ++d_smtCalls;
      if (res.isUnsat())
      {
        break;
      }
      if (!res.isSat())
      {
        throw std::runtime_error("solver returned unknown during enumeration");
      }

      ++total;

      std::vector<cvc5::Term> equalities;
      equalities.reserve(d_projectionVars.size());
      for (const cvc5::Term& var : d_projectionVars)
      {
        cvc5::Term value = d_solver.getValue(var);
        equalities.push_back(tm.mkTerm(cvc5::Kind::EQUAL, {var, value}));
      }

      cvc5::Term block;
      if (equalities.size() == 1)
      {
        block = tm.mkTerm(cvc5::Kind::NOT, {equalities[0]});
      }
      else
      {
        cvc5::Term conjunction = tm.mkTerm(cvc5::Kind::AND, equalities);
        block = tm.mkTerm(cvc5::Kind::NOT, {conjunction});
      }
      d_solver.assertFormula(block);
    }

    d_solver.pop();
    return total;
  }
  catch (...)
  {
    d_solver.pop();
    throw;
  }
}

long double ProjectedEnumerator::countWeighted(
    const std::unordered_map<cvc5::Term, TTCParser::LiteralWeight>& weights)
{
  d_smtCalls = 0;
  d_solver.push();
  long double total = 0.0L;
  try
  {
    if (d_projectionVars.empty())
    {
      cvc5::Result res = d_solver.checkSat();
      ++d_smtCalls;
      if (res.isSat())
      {
        total = 1.0L;
      }
      else if (res.isUnsat())
      {
        total = 0.0L;
      }
      else
      {
        throw std::runtime_error("solver returned unknown during enumeration");
      }
      d_solver.pop();
      return total;
    }

    auto& tm = ttc::getTermBuilder(d_solver);
    while (true)
    {
      cvc5::Result res = d_solver.checkSat();
      ++d_smtCalls;
      if (res.isUnsat())
      {
        break;
      }
      if (!res.isSat())
      {
        throw std::runtime_error("solver returned unknown during enumeration");
      }

      long double modelWeight = 1.0L;
      std::vector<cvc5::Term> equalities;
      equalities.reserve(d_projectionVars.size());
      for (const cvc5::Term& var : d_projectionVars)
      {
        cvc5::Term value = d_solver.getValue(var);
        equalities.push_back(tm.mkTerm(cvc5::Kind::EQUAL, {var, value}));
        if (var.getSort().isBoolean())
        {
          auto it = weights.find(var);
          const TTCParser::LiteralWeight defaultWeight;
          const TTCParser::LiteralWeight& weight =
              it == weights.end() ? defaultWeight : it->second;
          modelWeight *= value.getBooleanValue()
                             ? static_cast<long double>(weight.positive)
                             : static_cast<long double>(weight.negative);
        }
      }
      total += modelWeight;

      cvc5::Term block;
      if (equalities.size() == 1)
      {
        block = tm.mkTerm(cvc5::Kind::NOT, {equalities[0]});
      }
      else
      {
        cvc5::Term conjunction = tm.mkTerm(cvc5::Kind::AND, equalities);
        block = tm.mkTerm(cvc5::Kind::NOT, {conjunction});
      }
      d_solver.assertFormula(block);
    }

    d_solver.pop();
    return total;
  }
  catch (...)
  {
    d_solver.pop();
    throw;
  }
}
