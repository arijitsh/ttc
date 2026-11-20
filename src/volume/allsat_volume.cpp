#include "volume/allsat_volume.hpp"

#include <algorithm>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <cstdlib>
#include <gmpxx.h>

#include "Aiger/AigerParser.hpp"
#include "AllSatAlgo/Blocking/TseitinEnc/AllSatAlgoTseitinEnc.hpp"
#include "Globals/AllSatSolverGloblas.hpp"
#include "Globals/TernaryVal.hpp"
#include "Utilities/InputParser.hpp"

namespace ttc
{
namespace
{

class MutableAigerParser : public AigerParser
{
 public:
  void load(const BooleanAbstractionAigData& data)
  {
    m_Inputs.assign(data.inputs.begin(), data.inputs.end());
    m_Outputs.assign(data.outputs.begin(), data.outputs.end());
    m_AndGates.clear();
    m_AndGates.reserve(data.andGates.size());
    m_IsVarRef.assign(static_cast<std::size_t>(data.maxVariable) + 1, false);
    if (!m_IsVarRef.empty())
    {
      m_IsVarRef[0] = true;
    }
    auto markVar = [&](std::uint32_t lit) {
      std::uint32_t var = lit >> 1U;
      if (var >= m_IsVarRef.size())
      {
        m_IsVarRef.resize(var + 1, false);
      }
      m_IsVarRef[var] = true;
    };
    for (std::uint32_t lit : data.inputs)
    {
      markVar(lit);
    }
    for (const auto& gate : data.andGates)
    {
      markVar(gate.lhs);
      markVar(gate.rhs0);
      markVar(gate.rhs1);
      m_AndGates.emplace_back(gate.lhs, gate.rhs0, gate.rhs1);
    }
  }
};

class AllSatDnfEnumerator : public AllSatAlgoTseitinEnc
{
 public:
  AllSatDnfEnumerator(const BooleanAbstractionAigData& data, const InputParser& parser)
      : AllSatAlgoTseitinEnc(parser)
  {
    initialize(data);
  }

  std::vector<INPUT_ASSIGNMENT> enumerate()
  {
    std::vector<INPUT_ASSIGNMENT> assignments;
    SOLVER_RET_STATUS status = m_Solver->Solve();
    while (status == SAT_RET_STATUS)
    {
      INPUT_ASSIGNMENT model = m_Solver->GetAssignmentForAIGLits(m_Inputs);
      INPUT_ASSIGNMENT generalized = GeneralizeModel(model);
      if (m_IsTimeOut)
      {
        break;
      }
      assignments.push_back(generalized);
      BlockModel(generalized);
      status = m_Solver->Solve();
    }
    return assignments;
  }

 private:
  void initialize(const BooleanAbstractionAigData& data)
  {
    MutableAigerParser parser;
    parser.load(data);
    m_AigParser = parser;
    if (m_CirSimulation != nullptr)
    {
      delete m_CirSimulation;
      m_CirSimulation = nullptr;
    }
    if (m_UseCirSim)
    {
      m_CirSimulation = new CirSim(m_AigParser,
                                   m_UseTopToBotSim ? SimStrat::TopToBot : SimStrat::BotToTop);
    }
    m_Inputs = m_AigParser.GetInputs();
    m_InputSize = m_Inputs.size();
    m_Solver->InitializeSolver(m_AigParser);
    if (m_UseDualSolver && m_DualSolver != nullptr)
    {
      m_DualSolver->InitializeSolver(m_AigParser);
    }
  }
};

struct LiteralPairHash
{
  std::size_t operator()(const std::pair<std::uint32_t, std::uint32_t>& value) const noexcept
  {
    return (static_cast<std::size_t>(value.first) << 32) ^ static_cast<std::size_t>(value.second);
  }
};

class ManualAigBuilder
{
 public:
  ManualAigBuilder() : d_nextLiteral(2) { d_data.maxVariable = 1; }

  BooleanAbstractionAigData build(const cvc5::Term& term)
  {
    d_data.inputs.clear();
    d_data.outputs.clear();
    d_data.andGates.clear();
    d_data.mapping.clear();
    d_data.maxVariable = 1;
    d_cache.clear();
    d_atomLiterals.clear();
    d_andCache.clear();
    d_nextLiteral = 2;

    std::uint32_t literal = visit(term);
    d_data.outputs.push_back(literal);

    std::sort(d_data.mapping.begin(), d_data.mapping.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    std::sort(d_data.inputs.begin(), d_data.inputs.end());

    return d_data;
  }

 private:
  static std::uint32_t negate(std::uint32_t lit) { return lit ^ 1U; }

  std::uint32_t allocateLiteral()
  {
    std::uint32_t literal = d_nextLiteral;
    d_nextLiteral += 2;
    d_data.maxVariable = std::max(d_data.maxVariable, literal >> 1U);
    return literal;
  }

  std::uint32_t getOrCreateAtom(const cvc5::Term& term)
  {
    auto it = d_atomLiterals.find(term);
    if (it != d_atomLiterals.end())
    {
      return it->second;
    }
    std::uint32_t literal = allocateLiteral();
    d_atomLiterals.emplace(term, literal);
    d_data.inputs.push_back(literal);
    d_data.mapping.emplace_back(literal >> 1U, term);
    return literal;
  }

  std::uint32_t makeAnd(std::uint32_t lhs, std::uint32_t rhs)
  {
    if (lhs == 0U || rhs == 0U)
    {
      return 0U;
    }
    if (lhs == 1U)
    {
      return rhs;
    }
    if (rhs == 1U)
    {
      return lhs;
    }
    if (lhs == rhs)
    {
      return lhs;
    }
    if (lhs == negate(rhs))
    {
      return 0U;
    }
    if (lhs > rhs)
    {
      std::swap(lhs, rhs);
    }
    auto key = std::make_pair(lhs, rhs);
    auto it = d_andCache.find(key);
    if (it != d_andCache.end())
    {
      return it->second;
    }
    std::uint32_t literal = allocateLiteral();
    d_data.andGates.push_back({literal, lhs, rhs});
    d_andCache.emplace(key, literal);
    return literal;
  }

  std::uint32_t makeOr(std::uint32_t lhs, std::uint32_t rhs)
  {
    if (lhs == 1U || rhs == 1U)
    {
      return 1U;
    }
    if (lhs == 0U)
    {
      return rhs;
    }
    if (rhs == 0U)
    {
      return lhs;
    }
    if (lhs == rhs)
    {
      return lhs;
    }
    if (lhs == negate(rhs))
    {
      return 1U;
    }
    std::uint32_t andLit = makeAnd(negate(lhs), negate(rhs));
    return negate(andLit);
  }

  std::uint32_t visit(const cvc5::Term& term)
  {
    auto cached = d_cache.find(term);
    if (cached != d_cache.end())
    {
      return cached->second;
    }

    cvc5::Kind kind = term.getKind();
    std::uint32_t result = 0U;
    switch (kind)
    {
      case cvc5::Kind::CONST_BOOLEAN:
        result = term.getBooleanValue() ? 1U : 0U;
        break;
      case cvc5::Kind::NOT:
        result = negate(visit(term[0]));
        break;
      case cvc5::Kind::AND:
      {
        if (term.getNumChildren() == 0)
        {
          result = 1U;
        }
        else
        {
          result = visit(term[0]);
          for (std::size_t i = 1; i < term.getNumChildren(); ++i)
          {
            result = makeAnd(result, visit(term[i]));
          }
        }
        break;
      }
      case cvc5::Kind::OR:
      {
        if (term.getNumChildren() == 0)
        {
          result = 0U;
        }
        else
        {
          result = visit(term[0]);
          for (std::size_t i = 1; i < term.getNumChildren(); ++i)
          {
            result = makeOr(result, visit(term[i]));
          }
        }
        break;
      }
      case cvc5::Kind::IMPLIES:
      {
        std::uint32_t lhs = visit(term[0]);
        std::uint32_t rhs = visit(term[1]);
        result = makeOr(negate(lhs), rhs);
        break;
      }
      case cvc5::Kind::ITE:
      {
        std::uint32_t cond = visit(term[0]);
        std::uint32_t thenLit = visit(term[1]);
        std::uint32_t elseLit = visit(term[2]);
        std::uint32_t positive = makeAnd(cond, thenLit);
        std::uint32_t negative = makeAnd(negate(cond), elseLit);
        result = makeOr(positive, negative);
        break;
      }
      case cvc5::Kind::XOR:
      {
        std::uint32_t a = visit(term[0]);
        std::uint32_t b = visit(term[1]);
        std::uint32_t left = makeAnd(negate(a), b);
        std::uint32_t right = makeAnd(a, negate(b));
        result = makeOr(left, right);
        break;
      }
      case cvc5::Kind::EQUAL:
      {
        if (term[0].getSort().isBoolean())
        {
          std::uint32_t a = visit(term[0]);
          std::uint32_t b = visit(term[1]);
          std::uint32_t pos = makeAnd(a, b);
          std::uint32_t neg = makeAnd(negate(a), negate(b));
          result = makeOr(pos, neg);
        }
        else
        {
          result = getOrCreateAtom(term);
        }
        break;
      }
      default:
      {
        if (!term.getSort().isBoolean())
        {
          result = getOrCreateAtom(term);
        }
        else
        {
          result = getOrCreateAtom(term);
        }
        break;
      }
    }

    if (kind != cvc5::Kind::CONST_BOOLEAN && term.getSort().isBoolean())
    {
      d_cache.emplace(term, result);
    }
    return result;
  }

  BooleanAbstractionAigData d_data;
  std::unordered_map<cvc5::Term, std::uint32_t> d_cache;
  std::unordered_map<cvc5::Term, std::uint32_t> d_atomLiterals;
  std::unordered_map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t, LiteralPairHash> d_andCache;
  std::uint32_t d_nextLiteral;
};


}  // namespace

AllSatVolumeResult enumeratePolytopes(const BooleanAbstractionAigData& data,
                                      cvc5::Solver& solver)
{
  AllSatVolumeResult result;
  if (data.outputs.empty())
  {
    return result;
  }

  int argc = 1;
  char arg0[] = "allsat";
  char* argv[] = {arg0};
  InputParser parser(argc, argv);
  parser.AppendParams({"/general/print_enumer", "0",
                       "/mode", "roc",
                       "/alg/blocking/use_ucore", "1",
                       "/alg/blocking/use_ipasir_for_dual", "0"
                      });

  std::cout << "c starting DNFization" << std::endl;
  AllSatDnfEnumerator enumerator(data, parser);
  std::vector<INPUT_ASSIGNMENT> rawAssignments = enumerator.enumerate();
  std::cout << "c done DNFization" << std::endl;

  std::unordered_map<std::uint32_t, cvc5::Term> varToTerm;
  varToTerm.reserve(data.mapping.size());
  for (const auto& entry : data.mapping)
  {
    varToTerm.emplace(entry.first, entry.second);
  }

  auto& tm = ttc::getTermBuilder(solver);
  std::unordered_map<std::string, std::pair<Polytope, std::size_t>> uniquePolytopes;
  std::vector<std::string> keyOrder;
  keyOrder.reserve(rawAssignments.size());

  for (const auto& assignment : rawAssignments)
  {
    Polytope poly;
    poly.booleanLiterals.reserve(assignment.size());
    std::unordered_map<cvc5::Term, std::optional<bool>> termAssignments;
    termAssignments.reserve(assignment.size());
    bool contradictory = false;

    for (const auto& entry : assignment)
    {
      std::uint32_t lit = entry.first;
      TVal value = entry.second;
      std::uint32_t var = lit >> 1U;

      std::optional<bool> literalValue;
      if (value == TVal::True)
      {
        literalValue = true;
      }
      else if (value == TVal::False)
      {
        literalValue = false;
      }
      else if (value == TVal::DontCare)
      {
        literalValue = std::nullopt;
      }
      else
      {
        continue;
      }

      poly.booleanLiterals.push_back({var, literalValue});

      auto termIt = varToTerm.find(var);
      if (termIt != varToTerm.end())
      {
        auto [it, inserted] = termAssignments.emplace(termIt->second, literalValue);
        if (!inserted)
        {
          auto& existing = it->second;
          if (literalValue.has_value())
          {
            if (!existing.has_value())
            {
              existing = literalValue;
            }
            else if (existing.value() != literalValue.value())
            {
              contradictory = true;
              break;
            }
          }
        }
      }
    }

    if (contradictory)
    {
      continue;
    }

    std::vector<cvc5::Term> conjuncts;
    std::vector<std::string> keyParts;
    conjuncts.reserve(termAssignments.size());
    poly.termLiterals.reserve(termAssignments.size());
    std::size_t specifiedCount = 0;

    for (const auto& entry : termAssignments)
    {
      const cvc5::Term& term = entry.first;
      const std::optional<bool>& value = entry.second;
      poly.termLiterals.push_back({term, value});
      if (value.has_value())
      {
        ++specifiedCount;
        if (value.value())
        {
          conjuncts.push_back(term);
          keyParts.push_back(term.toString() + "=1");
        }
        else
        {
          conjuncts.push_back(tm.mkTerm(cvc5::Kind::NOT, {term}));
          keyParts.push_back(term.toString() + "=0");
        }
      }
    }

    std::sort(keyParts.begin(), keyParts.end());
    std::ostringstream keyBuilder;
    for (const auto& part : keyParts)
    {
      keyBuilder << part << ';';
    }
    std::string key = keyBuilder.str();

    if (conjuncts.empty())
    {
      poly.formula = tm.mkBoolean(true);
    }
    else if (conjuncts.size() == 1)
    {
      poly.formula = conjuncts.front();
    }
    else
    {
      poly.formula = tm.mkTerm(cvc5::Kind::AND, conjuncts);
    }

    auto it = uniquePolytopes.find(key);
    if (it == uniquePolytopes.end())
    {
      keyOrder.push_back(key);
      uniquePolytopes.emplace(key, std::make_pair(poly, specifiedCount));
    }
    else if (specifiedCount < it->second.second)
    {
      it->second.first = poly;
      it->second.second = specifiedCount;
    }
  }

  std::vector<std::vector<std::string>> keyTerms;
  keyTerms.reserve(keyOrder.size());
  for (const auto& key : keyOrder)
  {
    std::vector<std::string> terms;
    std::stringstream ss(key);
    std::string token;
    while (std::getline(ss, token, ';'))
    {
      if (!token.empty())
      {
        terms.push_back(token);
      }
    }
    keyTerms.push_back(std::move(terms));
  }

  auto isSubset = [](const std::vector<std::string>& smaller,
                     const std::vector<std::string>& larger) {
    if (smaller.size() > larger.size())
    {
      return false;
    }
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < smaller.size() && j < larger.size())
    {
      if (smaller[i] == larger[j])
      {
        ++i;
        ++j;
      }
      else if (smaller[i] > larger[j])
      {
        ++j;
      }
      else
      {
        return false;
      }
    }
    return i == smaller.size();
  };

  std::vector<bool> keep(keyOrder.size(), true);
  for (std::size_t i = 0; i < keyOrder.size(); ++i)
  {
    if (!keep[i])
    {
      continue;
    }
    for (std::size_t j = 0; j < keyOrder.size(); ++j)
    {
      if (i == j || !keep[j])
      {
        continue;
      }
      if (keyTerms[j].size() <= keyTerms[i].size() &&
          isSubset(keyTerms[j], keyTerms[i]) &&
          keyTerms[j].size() < keyTerms[i].size())
      {
        keep[i] = false;
        break;
      }
    }
  }

  std::vector<Polytope> candidates;
  candidates.reserve(uniquePolytopes.size());
  for (std::size_t idx = 0; idx < keyOrder.size(); ++idx)
  {
    if (!keep[idx])
    {
      continue;
    }
    const auto& key = keyOrder[idx];
    auto it = uniquePolytopes.find(key);
    if (it != uniquePolytopes.end())
    {
      candidates.push_back(std::move(it->second.first));
    }
  }

  result.polytopes = std::move(candidates);
  return result;
}

BooleanAbstractionAigData buildBooleanAigFromTerm(const cvc5::Term& term)
{
  ManualAigBuilder builder;
  return builder.build(term);
}

}  // namespace ttc

