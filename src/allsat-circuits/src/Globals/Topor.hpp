#pragma once

#include <cstdint>
#include <cstdlib>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cvc5/cvc5.h>

#include "Globals/AllSatGloblas.hpp"

namespace Topor
{

enum class TToporReturnVal
{
  RET_SAT,
  RET_UNSAT,
  RET_TIMEOUT_LOCAL,
  RET_CONFLICT_OUT,
  RET_MEM_OUT,
  RET_USER_INTERRUPT,
  RET_INDEX_TOO_NARROW,
  RET_PARAM_ERROR,
  RET_ASSUMPTION_REQUIRED_ERROR,
  RET_TIMEOUT_GLOBAL,
  RET_DRAT_FILE_PROBLEM,
  RET_EXOTIC_ERROR
};

enum class TToporLitVal
{
  VAL_UNDEF,
  VAL_SATISFIED,
  VAL_UNSATISFIED
};

template <typename LitT, typename IndexT, bool Compress>
class CTopor
{
 public:
  CTopor()
      : d_solver(),
        d_lastResult(TToporReturnVal::RET_SAT)
  {
    d_solver.setOption("produce-models", "true");
    d_solver.setOption("produce-unsat-cores", "true");
    d_solver.setOption("incremental", "true");
    d_varTerms.resize(2);
  }

  void SetParam(const std::string& name, double value)
  {
    if (name == "/timeout/global")
    {
      d_timeoutSeconds = value;
    }
    (void)value;
  }

  void AddClause(const std::vector<LitT>& clause)
  {
    bool clauseTrue = false;
    std::vector<cvc5::Term> disj;
    disj.reserve(clause.size());
    for (LitT lit : clause)
    {
      if (lit == CONST_LIT_TRUE)
      {
        clauseTrue = true;
        break;
      }
      if (lit == CONST_LIT_FALSE)
      {
        continue;
      }
      std::uint32_t var = static_cast<std::uint32_t>(std::abs(lit));
      ensureVar(var);
      cvc5::Term term = d_varTerms[var];
      if (lit < 0)
      {
        term = d_solver.mkTerm(cvc5::Kind::NOT, {term});
      }
      disj.push_back(term);
    }
    if (clauseTrue)
    {
      return;
    }
    if (disj.empty())
    {
      d_solver.assertFormula(d_solver.mkBoolean(false));
    }
    else if (disj.size() == 1)
    {
      d_solver.assertFormula(disj.front());
    }
    else
    {
      d_solver.assertFormula(d_solver.mkTerm(cvc5::Kind::OR, disj));
    }
  }

  TToporReturnVal Solve()
  {
    return solveInternal({});
  }

  TToporReturnVal Solve(const std::vector<LitT>& assumptions)
  {
    return solveInternal(assumptions);
  }

  void FixPolarity(LitT) {}

  void BoostScore(LitT) {}

  TToporLitVal GetLitValue(LitT lit) const
  {
    if (lit == CONST_LIT_TRUE)
    {
      return TToporLitVal::VAL_SATISFIED;
    }
    if (lit == CONST_LIT_FALSE)
    {
      return TToporLitVal::VAL_UNSATISFIED;
    }
    std::uint32_t var = static_cast<std::uint32_t>(std::abs(lit));
    auto it = d_lastAssignment.find(var);
    if (it == d_lastAssignment.end())
    {
      return TToporLitVal::VAL_UNDEF;
    }
    bool value = it->second;
    bool satisfied = (lit > 0) ? value : !value;
    return satisfied ? TToporLitVal::VAL_SATISFIED : TToporLitVal::VAL_UNSATISFIED;
  }

  bool IsAssumptionRequired(std::size_t pos) const
  {
    if (pos >= d_assumptionRequired.size())
    {
      return false;
    }
    return d_assumptionRequired[pos];
  }

  void SetConflictLimit(int) {}

 private:
  void ensureVar(std::uint32_t var)
  {
    if (var >= d_varTerms.size())
    {
      d_varTerms.resize(var + 1);
    }
    if (d_varTerms[var].isNull())
    {
      std::ostringstream name;
      name << "v" << var;
      d_varTerms[var] = d_solver.mkConst(d_solver.getBooleanSort(), name.str());
    }
  }

  void refreshAssignment()
  {
    d_lastAssignment.clear();
    for (std::size_t idx = 1; idx < d_varTerms.size(); ++idx)
    {
      const cvc5::Term& term = d_varTerms[idx];
      if (term.isNull())
      {
        continue;
      }
      cvc5::Term value = d_solver.getValue(term);
      d_lastAssignment[static_cast<std::uint32_t>(idx)] = value.getBooleanValue();
    }
  }

  TToporReturnVal solveInternal(const std::vector<LitT>& assumptions)
  {
    std::vector<cvc5::Term> assumptionTerms;
    assumptionTerms.reserve(assumptions.size());
    std::vector<std::size_t> indices;
    d_assumptionRequired.assign(assumptions.size(), false);

    for (std::size_t i = 0; i < assumptions.size(); ++i)
    {
      LitT lit = assumptions[i];
      if (lit == CONST_LIT_TRUE)
      {
        continue;
      }
      if (lit == CONST_LIT_FALSE)
      {
        d_assumptionRequired[i] = true;
        d_lastResult = TToporReturnVal::RET_UNSAT;
        return d_lastResult;
      }
      std::uint32_t var = static_cast<std::uint32_t>(std::abs(lit));
      ensureVar(var);
      cvc5::Term term = d_varTerms[var];
      if (lit < 0)
      {
        term = d_solver.mkTerm(cvc5::Kind::NOT, {term});
      }
      indices.push_back(i);
      assumptionTerms.push_back(term);
    }

    cvc5::Result result;
    if (assumptionTerms.empty())
    {
      result = d_solver.checkSat();
    }
    else
    {
      result = d_solver.checkSatAssuming(assumptionTerms);
    }

    if (result.isSat())
    {
      refreshAssignment();
      d_lastResult = TToporReturnVal::RET_SAT;
      return d_lastResult;
    }
    if (result.isUnsat())
    {
      auto core = d_solver.getUnsatCore();
      for (const cvc5::Term& term : core)
      {
        for (std::size_t j = 0; j < assumptionTerms.size(); ++j)
        {
          if (term == assumptionTerms[j])
          {
            d_assumptionRequired[indices[j]] = true;
          }
        }
      }
      d_lastResult = TToporReturnVal::RET_UNSAT;
      return d_lastResult;
    }

    d_lastResult = TToporReturnVal::RET_TIMEOUT_GLOBAL;
    return d_lastResult;
  }

  cvc5::Solver d_solver;
  std::vector<cvc5::Term> d_varTerms;
  std::unordered_map<std::uint32_t, bool> d_lastAssignment;
  std::vector<bool> d_assumptionRequired;
  TToporReturnVal d_lastResult;
  double d_timeoutSeconds = -1.0;
};

}  // namespace Topor

