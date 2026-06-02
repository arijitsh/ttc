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
#include <isl/constraint.h>
#include <isl/ctx.h>
#include <isl/local_space.h>
#include <isl/set.h>
#include <isl/space.h>
#include <isl/val.h>

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

using TermIndexMap = std::unordered_map<cvc5::Term, std::size_t>;

class SolverPushGuard
{
 public:
  explicit SolverPushGuard(cvc5::Solver& solver) : d_solver(solver), d_active(true)
  {
    d_solver.push();
  }

  ~SolverPushGuard()
  {
    if (d_active)
    {
      d_solver.pop();
    }
  }

  SolverPushGuard(const SolverPushGuard&) = delete;
  SolverPushGuard& operator=(const SolverPushGuard&) = delete;

  void dismiss()
  {
    if (d_active)
    {
      d_solver.pop();
      d_active = false;
    }
  }

 private:
  cvc5::Solver& d_solver;
  bool d_active;
};

mpq_class parseRationalString(const std::string& value)
{
  if (value.empty())
  {
    return mpq_class(0);
  }
  return mpq_class(value);
}

std::optional<mpq_class> getConstantValue(const cvc5::Term& term)
{
  if (term.isRealValue())
  {
    return parseRationalString(term.getRealValue());
  }
  if (term.isIntegerValue())
  {
    return parseRationalString(term.getIntegerValue());
  }
  switch (term.getKind())
  {
    case cvc5::Kind::TO_REAL:
    case cvc5::Kind::TO_INTEGER:
      return getConstantValue(term[0]);
    case cvc5::Kind::NEG:
    {
      if (auto inner = getConstantValue(term[0]))
      {
        return -(*inner);
      }
      return std::nullopt;
    }
    case cvc5::Kind::SUB:
      if (term.getNumChildren() == 1)
      {
        if (auto inner = getConstantValue(term[0]))
        {
          return -(*inner);
        }
        return std::nullopt;
      }
      break;
    case cvc5::Kind::DIVISION:
      if (term.getNumChildren() == 2)
      {
        auto numerator = getConstantValue(term[0]);
        auto denominator = getConstantValue(term[1]);
        if (numerator && denominator && *denominator != 0)
        {
          return *numerator / *denominator;
        }
        return std::nullopt;
      }
      break;
    case cvc5::Kind::MULT:
    {
      mpq_class product(1);
      bool hasFactor = false;
      for (const auto& child : term)
      {
        auto factor = getConstantValue(child);
        if (!factor)
        {
          return std::nullopt;
        }
        product *= *factor;
        hasFactor = true;
      }
      if (hasFactor)
      {
        return product;
      }
      return std::nullopt;
    }
    default:
      break;
  }
  return std::nullopt;
}

bool accumulateLinear(const cvc5::Term& term,
                      const mpq_class& multiplier,
                      const TermIndexMap& index,
                      std::vector<mpq_class>& coefficients,
                      mpq_class& constant)
{
  if (multiplier == 0)
  {
    return true;
  }

  if (term.isRealValue())
  {
    constant += multiplier * parseRationalString(term.getRealValue());
    return true;
  }
  if (term.isIntegerValue())
  {
    constant += multiplier * parseRationalString(term.getIntegerValue());
    return true;
  }

  if (term.getNumChildren() == 0)
  {
    if (term.getSort().isReal() || term.getSort().isInteger())
    {
      auto it = index.find(term);
      if (it == index.end())
      {
        return false;
      }
      coefficients[it->second] += multiplier;
      return true;
    }
    return false;
  }

  switch (term.getKind())
  {
    case cvc5::Kind::ADD:
      for (const auto& child : term)
      {
        if (!accumulateLinear(child, multiplier, index, coefficients, constant))
        {
          return false;
        }
      }
      return true;
    case cvc5::Kind::SUB:
      if (term.getNumChildren() == 1)
      {
        return accumulateLinear(term[0], -multiplier, index, coefficients, constant);
      }
      if (!accumulateLinear(term[0], multiplier, index, coefficients, constant))
      {
        return false;
      }
      for (std::size_t i = 1, n = term.getNumChildren(); i < n; ++i)
      {
        if (!accumulateLinear(term[i], -multiplier, index, coefficients, constant))
        {
          return false;
        }
      }
      return true;
    case cvc5::Kind::NEG:
      return accumulateLinear(term[0], -multiplier, index, coefficients, constant);
    case cvc5::Kind::TO_REAL:
    case cvc5::Kind::TO_INTEGER:
      return accumulateLinear(term[0], multiplier, index, coefficients, constant);
    case cvc5::Kind::DIVISION:
      if (term.getNumChildren() == 2)
      {
        auto denom = getConstantValue(term[1]);
        if (!denom || *denom == 0)
        {
          return false;
        }
        return accumulateLinear(term[0], multiplier / *denom, index, coefficients,
                                constant);
      }
      return false;
    case cvc5::Kind::MULT:
    {
      mpq_class factor(1);
      std::optional<cvc5::Term> expr;
      for (const auto& child : term)
      {
        auto constantFactor = getConstantValue(child);
        if (constantFactor)
        {
          factor *= *constantFactor;
          continue;
        }
        if (expr.has_value())
        {
          return false;
        }
        expr = child;
      }
      if (!expr.has_value())
      {
        constant += multiplier * factor;
        return true;
      }
      return accumulateLinear(expr.value(), multiplier * factor, index, coefficients,
                              constant);
    }
    default:
      break;
  }
  return false;
}

struct LinearInequality
{
  std::vector<mpq_class> coefficients;
  mpq_class bound{0};
};

bool appendLessEqualConstraint(const cvc5::Term& lhs,
                               const cvc5::Term& rhs,
                               const TermIndexMap& index,
                               std::size_t dimension,
                               std::vector<LinearInequality>& inequalities,
                               bool& empty)
{
  std::vector<mpq_class> coeffs(dimension, mpq_class(0));
  mpq_class constant(0);
  if (!accumulateLinear(lhs, mpq_class(1), index, coeffs, constant) ||
      !accumulateLinear(rhs, mpq_class(-1), index, coeffs, constant))
  {
    return false;
  }
  bool allZero = std::all_of(coeffs.begin(), coeffs.end(),
                             [](const mpq_class& value) { return value == 0; });
  if (allZero)
  {
    if (constant <= 0)
    {
      return true;
    }
    empty = true;
    return true;
  }
  LinearInequality inequality;
  inequality.coefficients = std::move(coeffs);
  inequality.bound = -constant;
  inequalities.push_back(std::move(inequality));
  return true;
}

bool appendComparisonConstraints(const cvc5::Term& atom,
                                 bool positive,
                                 const TermIndexMap& index,
                                 std::size_t dimension,
                                 std::vector<LinearInequality>& inequalities,
                                 bool& empty)
{
  cvc5::Kind kind = atom.getKind();
  if (!positive)
  {
    switch (kind)
    {
      case cvc5::Kind::LEQ:
        kind = cvc5::Kind::GEQ;
        break;
      case cvc5::Kind::LT:
        kind = cvc5::Kind::GT;
        break;
      case cvc5::Kind::GEQ:
        kind = cvc5::Kind::LEQ;
        break;
      case cvc5::Kind::GT:
        kind = cvc5::Kind::LT;
        break;
      case cvc5::Kind::EQUAL:
        throw std::runtime_error(
            "negated equality is not supported during polytope simplification");
      default:
        break;
    }
  }

  switch (kind)
  {
    case cvc5::Kind::LEQ:
    case cvc5::Kind::LT:
      return appendLessEqualConstraint(atom[0], atom[1], index, dimension,
                                       inequalities, empty);
    case cvc5::Kind::GEQ:
    case cvc5::Kind::GT:
      return appendLessEqualConstraint(atom[1], atom[0], index, dimension,
                                       inequalities, empty);
    case cvc5::Kind::EQUAL:
      if (!appendLessEqualConstraint(atom[0], atom[1], index, dimension,
                                     inequalities, empty))
      {
        return false;
      }
      return appendLessEqualConstraint(atom[1], atom[0], index, dimension,
                                       inequalities, empty);
    default:
      break;
  }
  return false;
}

bool collectConstraints(const cvc5::Term& term,
                        bool positive,
                        const TermIndexMap& index,
                        std::size_t dimension,
                        std::vector<LinearInequality>& inequalities,
                        bool& empty)
{
  if (empty)
  {
    return true;
  }

  if (term.getKind() == cvc5::Kind::NOT)
  {
    return collectConstraints(term[0], !positive, index, dimension, inequalities,
                              empty);
  }

  if (term.isBooleanValue())
  {
    bool value = term.getBooleanValue();
    bool satisfied = positive ? value : !value;
    if (!satisfied)
    {
      empty = true;
    }
    return true;
  }

  if (term.getKind() == cvc5::Kind::AND)
  {
    for (const auto& child : term)
    {
      if (!collectConstraints(child, positive, index, dimension, inequalities, empty))
      {
        return false;
      }
    }
    return true;
  }

  return appendComparisonConstraints(term, positive, index, dimension, inequalities,
                                     empty);
}

void collectLinearVariables(const cvc5::Term& term,
                            std::unordered_set<cvc5::Term>& variables)
{
  if (term.isRealValue() || term.isIntegerValue())
  {
    return;
  }

  if (term.getNumChildren() == 0)
  {
    if (term.getSort().isReal() || term.getSort().isInteger())
    {
      variables.insert(term);
    }
    return;
  }

  for (const auto& child : term)
  {
    collectLinearVariables(child, variables);
  }
}

mpz_class gcdAbs(const mpz_class& a, const mpz_class& b)
{
  mpz_class result;
  mpz_gcd(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  if (result < 0)
  {
    result = -result;
  }
  return result;
}

mpz_class lcmAbs(const mpz_class& a, const mpz_class& b)
{
  if (a == 0)
  {
    return b < 0 ? -b : b;
  }
  if (b == 0)
  {
    return a < 0 ? -a : a;
  }
  mpz_class absA = a < 0 ? -a : a;
  mpz_class absB = b < 0 ? -b : b;
  mpz_class g = gcdAbs(absA, absB);
  return (absA / g) * absB;
}

mpz_class computeScale(const LinearInequality& inequality)
{
  mpz_class scale(1);
  for (const auto& coeff : inequality.coefficients)
  {
    mpz_class denom = coeff.get_den();
    if (denom < 0)
    {
      denom = -denom;
    }
    scale = lcmAbs(scale, denom);
  }
  mpz_class boundDenom = inequality.bound.get_den();
  if (boundDenom < 0)
  {
    boundDenom = -boundDenom;
  }
  scale = lcmAbs(scale, boundDenom);
  if (scale == 0)
  {
    scale = 1;
  }
  return scale;
}

struct IslBasicSetDeleter
{
  void operator()(isl_basic_set* set) const noexcept
  {
    isl_basic_set_free(set);
  }
};

std::unique_ptr<isl_basic_set, IslBasicSetDeleter> buildIslBasicSet(
    isl_ctx* ctx, const Polytope& poly)
{
  std::unordered_set<cvc5::Term> variables;
  for (const auto& literal : poly.termLiterals)
  {
    if (!literal.value.has_value())
    {
      continue;
    }
    collectLinearVariables(literal.term, variables);
  }

  std::vector<cvc5::Term> ordered;
  ordered.reserve(variables.size());
  for (const auto& term : variables)
  {
    ordered.push_back(term);
  }
  std::sort(ordered.begin(), ordered.end(),
            [](const cvc5::Term& a, const cvc5::Term& b) {
              return a.toString() < b.toString();
            });

  TermIndexMap index;
  index.reserve(ordered.size());
  for (std::size_t i = 0; i < ordered.size(); ++i)
  {
    index.emplace(ordered[i], i);
  }

  std::vector<LinearInequality> inequalities;
  inequalities.reserve(poly.termLiterals.size());
  bool empty = false;
  for (const auto& literal : poly.termLiterals)
  {
    if (!literal.value.has_value())
    {
      continue;
    }
    if (!collectConstraints(literal.term, literal.value.value(), index,
                            ordered.size(), inequalities, empty))
    {
      return nullptr;
    }
    if (empty)
    {
      break;
    }
  }

  if (empty)
  {
    return nullptr;
  }

  isl_space* space = isl_space_set_alloc(ctx, 0, ordered.size());
  isl_basic_set* set = isl_basic_set_universe(space);
  if (!set)
  {
    return nullptr;
  }

  if (!inequalities.empty())
  {
    isl_space* localSpace = isl_space_set_alloc(ctx, 0, ordered.size());
    isl_local_space* ls = isl_local_space_from_space(localSpace);
    if (!ls)
    {
      isl_basic_set_free(set);
      return nullptr;
    }

    for (const auto& inequality : inequalities)
    {
      mpz_class scale = computeScale(inequality);
      isl_constraint* constraint = isl_constraint_alloc_inequality(
          isl_local_space_copy(ls));
      if (!constraint)
      {
        isl_local_space_free(ls);
        isl_basic_set_free(set);
        return nullptr;
      }

      mpq_class scaledBound = inequality.bound * scale;
      mpz_class boundInt = scaledBound.get_num();
      std::string boundStr = boundInt.get_str();
      isl_val* boundVal = isl_val_read_from_str(ctx, boundStr.c_str());
      constraint = isl_constraint_set_constant_val(constraint, boundVal);

      for (std::size_t i = 0; i < inequality.coefficients.size(); ++i)
      {
        mpq_class scaledCoeff = inequality.coefficients[i] * scale;
        mpq_class negCoeff = -scaledCoeff;
        mpz_class coeffInt = negCoeff.get_num();
        std::string coeffStr = coeffInt.get_str();
        isl_val* coeffVal = isl_val_read_from_str(ctx, coeffStr.c_str());
        constraint = isl_constraint_set_coefficient_val(
            constraint, isl_dim_set, static_cast<int>(i), coeffVal);
      }

      set = isl_basic_set_add_constraint(set, constraint);
      if (!set)
      {
        isl_local_space_free(ls);
        return nullptr;
      }
    }
    isl_local_space_free(ls);
  }

  return std::unique_ptr<isl_basic_set, IslBasicSetDeleter>(set);
}

struct IslCtxDeleter
{
  void operator()(isl_ctx* ctx) const noexcept
  {
    if (ctx != nullptr)
    {
      isl_ctx_free(ctx);
    }
  }
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

  std::unique_ptr<isl_ctx, IslCtxDeleter> ctx(isl_ctx_alloc(), IslCtxDeleter{});
  result.polytopes.reserve(candidates.size());
  for (auto& poly : candidates)
  {
    if (!ctx)
    {
      result.polytopes.push_back(std::move(poly));
      continue;
    }

    auto islSet = buildIslBasicSet(ctx.get(), poly);
    if (!islSet)
    {
      result.polytopes.push_back(std::move(poly));
      continue;
    }

    isl_basic_set* raw = islSet.release();
    raw = isl_basic_set_detect_equalities(raw);
    raw = isl_basic_set_remove_redundancies(raw);
    islSet.reset(raw);
    if (!islSet)
    {
      result.polytopes.push_back(std::move(poly));
      continue;
    }

    isl_bool emptySet = isl_basic_set_is_empty(islSet.get());
    poly.islChecked = true;
    if (emptySet == isl_bool_true)
    {
      poly.islNonEmpty = false;
      try
      {
        if (!poly.formula.isNull())
        {
          SolverPushGuard guard(solver);
          solver.assertFormula(poly.formula);
          cvc5::Result res = solver.checkSat();
          if (res.isSat())
          {
            poly.islNonEmpty = true;
          }
        }
      }
      catch (const cvc5::CVC5ApiException&)
      {
        poly.islChecked = false;
        poly.islNonEmpty = true;
      }
    }
    else if (emptySet == isl_bool_false)
    {
      poly.islNonEmpty = true;
    }
    else
    {
      poly.islChecked = false;
      poly.islNonEmpty = true;
    }

    char* repr = isl_basic_set_to_str(islSet.get());
    if (repr != nullptr)
    {
      free(repr);
    }

    result.polytopes.push_back(std::move(poly));
  }

  return result;
}

BooleanAbstractionAigData buildBooleanAigFromTerm(const cvc5::Term& term)
{
  ManualAigBuilder builder;
  return builder.build(term);
}

}  // namespace ttc

