#include "parser.hpp"

#include <cctype>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "cache_mode.hpp"
#include "features.hpp"
#include "logger.hpp"
#include "stp_sat.hpp"
#include "volume/allsat_volume.hpp"

#include "arjun.hpp"
#ifdef TTC_ENABLE_DDNNF
#include "proj_ddnnf.hpp"
#include "var_order.hpp"
#endif

namespace {

struct ProjectionDeclarationParseResult
{
  std::string smtWithoutProjDecls;
  std::vector<std::string> projectionSymbols;
  struct WeightDeclaration
  {
    std::string symbol;
    bool negative = false;
    double weight = 0.0;
  };
  std::vector<WeightDeclaration> weightDeclarations;
  bool sawProjectionDeclaration = false;
  bool sawWeightDeclaration = false;
};

bool isSymbolDelimiter(char ch)
{
  return std::isspace(static_cast<unsigned char>(ch)) || ch == '(' || ch == ')';
}

std::size_t skipWhitespaceAndComments(const std::string& input,
                                      std::size_t pos,
                                      std::size_t end)
{
  while (pos < end)
  {
    if (std::isspace(static_cast<unsigned char>(input[pos])))
    {
      ++pos;
      continue;
    }
    if (input[pos] == ';')
    {
      while (pos < end && input[pos] != '\n')
      {
        ++pos;
      }
      continue;
    }
    break;
  }
  return pos;
}

std::string parseSmtSymbol(const std::string& input,
                           std::size_t& pos,
                           std::size_t end)
{
  if (pos >= end)
  {
    return {};
  }
  if (input[pos] == '|')
  {
    ++pos;
    std::string symbol;
    while (pos < end)
    {
      if (input[pos] == '\\' && pos + 1 < end)
      {
        symbol.push_back(input[pos + 1]);
        pos += 2;
        continue;
      }
      if (input[pos] == '|')
      {
        ++pos;
        return symbol;
      }
      symbol.push_back(input[pos++]);
    }
    throw std::runtime_error("unterminated quoted symbol in declare-projvar");
  }

  std::size_t start = pos;
  while (pos < end && !isSymbolDelimiter(input[pos]))
  {
    ++pos;
  }
  return input.substr(start, pos - start);
}

std::string parseSmtToken(const std::string& input,
                          std::size_t& pos,
                          std::size_t end)
{
  return parseSmtSymbol(input, pos, end);
}

std::size_t findMatchingCommandParen(const std::string& input,
                                     std::size_t openParen)
{
  std::size_t depth = 0;
  for (std::size_t i = openParen; i < input.size(); ++i)
  {
    char ch = input[i];
    if (ch == ';')
    {
      while (i < input.size() && input[i] != '\n')
      {
        ++i;
      }
      if (i >= input.size())
      {
        break;
      }
      continue;
    }
    if (ch == '|')
    {
      ++i;
      while (i < input.size())
      {
        if (input[i] == '\\' && i + 1 < input.size())
        {
          i += 2;
          continue;
        }
        if (input[i] == '|')
        {
          break;
        }
        ++i;
      }
      if (i >= input.size())
      {
        throw std::runtime_error("unterminated quoted symbol in SMT-LIB input");
      }
      continue;
    }
    if (ch == '"')
    {
      ++i;
      while (i < input.size())
      {
        if (input[i] == '"' && i + 1 < input.size() && input[i + 1] == '"')
        {
          i += 2;
          continue;
        }
        if (input[i] == '"')
        {
          break;
        }
        ++i;
      }
      if (i >= input.size())
      {
        throw std::runtime_error("unterminated string literal in SMT-LIB input");
      }
      continue;
    }
    if (ch == '(')
    {
      ++depth;
      continue;
    }
    if (ch == ')')
    {
      if (depth == 0)
      {
        throw std::runtime_error("unbalanced ')' in SMT-LIB input");
      }
      --depth;
      if (depth == 0)
      {
        return i + 1;
      }
    }
  }
  throw std::runtime_error("unterminated command in SMT-LIB input");
}

ProjectionDeclarationParseResult extractExtensionDeclarations(
    const std::string& smtFormula)
{
  ProjectionDeclarationParseResult result;
  result.smtWithoutProjDecls = smtFormula;
  std::size_t pos = 0;
  while (pos < smtFormula.size())
  {
    if (smtFormula[pos] == ';')
    {
      while (pos < smtFormula.size() && smtFormula[pos] != '\n')
      {
        ++pos;
      }
      continue;
    }
    if (smtFormula[pos] == '|')
    {
      ++pos;
      while (pos < smtFormula.size())
      {
        if (smtFormula[pos] == '\\' && pos + 1 < smtFormula.size())
        {
          pos += 2;
          continue;
        }
        if (smtFormula[pos++] == '|')
        {
          break;
        }
      }
      continue;
    }
    if (smtFormula[pos] != '(')
    {
      ++pos;
      continue;
    }

    std::size_t commandStart = pos;
    std::size_t commandEnd = findMatchingCommandParen(smtFormula, commandStart);
    std::size_t headPos =
        skipWhitespaceAndComments(smtFormula, commandStart + 1, commandEnd - 1);
    std::string head = parseSmtSymbol(smtFormula, headPos, commandEnd - 1);
    if (head == "declare-projvar")
    {
      result.sawProjectionDeclaration = true;
      while (true)
      {
        headPos = skipWhitespaceAndComments(smtFormula, headPos, commandEnd - 1);
        if (headPos >= commandEnd - 1 || smtFormula[headPos] == ')')
        {
          break;
        }
        if (smtFormula[headPos] == '(')
        {
          throw std::runtime_error(
              "declare-projvar expects variable symbols, not terms");
        }
        std::string symbol = parseSmtSymbol(smtFormula, headPos, commandEnd - 1);
        if (symbol.empty())
        {
          throw std::runtime_error("empty symbol in declare-projvar");
        }
        result.projectionSymbols.push_back(std::move(symbol));
      }
      std::fill(result.smtWithoutProjDecls.begin() + commandStart,
                result.smtWithoutProjDecls.begin() + commandEnd,
                ' ');
    }
    else if (head == "declare-weight")
    {
      result.sawWeightDeclaration = true;
      headPos = skipWhitespaceAndComments(smtFormula, headPos, commandEnd - 1);
      if (headPos >= commandEnd - 1 || smtFormula[headPos] == ')')
      {
        throw std::runtime_error("declare-weight expects a literal and weight");
      }

      ProjectionDeclarationParseResult::WeightDeclaration decl;
      if (smtFormula[headPos] == '(')
      {
        std::size_t litEnd = findMatchingCommandParen(smtFormula, headPos);
        if (litEnd > commandEnd - 1)
        {
          throw std::runtime_error("malformed literal in declare-weight");
        }
        std::size_t litPos =
            skipWhitespaceAndComments(smtFormula, headPos + 1, litEnd - 1);
        std::string litHead = parseSmtSymbol(smtFormula, litPos, litEnd - 1);
        if (litHead != "not" && litHead != "-")
        {
          throw std::runtime_error(
              "declare-weight supports only Boolean literals");
        }
        decl.negative = true;
        litPos = skipWhitespaceAndComments(smtFormula, litPos, litEnd - 1);
        decl.symbol = parseSmtSymbol(smtFormula, litPos, litEnd - 1);
        litPos = skipWhitespaceAndComments(smtFormula, litPos, litEnd - 1);
        if (decl.symbol.empty() || litPos < litEnd - 1)
        {
          throw std::runtime_error("malformed literal in declare-weight");
        }
        headPos = litEnd;
      }
      else
      {
        std::string literal = parseSmtSymbol(smtFormula, headPos, commandEnd - 1);
        if (!literal.empty() && literal[0] == '-' && literal.size() > 1)
        {
          decl.negative = true;
          decl.symbol = literal.substr(1);
        }
        else
        {
          decl.symbol = std::move(literal);
        }
      }

      headPos = skipWhitespaceAndComments(smtFormula, headPos, commandEnd - 1);
      std::string weightToken = parseSmtToken(smtFormula, headPos, commandEnd - 1);
      if (decl.symbol.empty() || weightToken.empty())
      {
        throw std::runtime_error("declare-weight expects a literal and weight");
      }
      std::size_t parsedChars = 0;
      try
      {
        decl.weight = std::stod(weightToken, &parsedChars);
      }
      catch (const std::exception&)
      {
        throw std::runtime_error("invalid weight in declare-weight: '" +
                                 weightToken + "'");
      }
      if (parsedChars != weightToken.size())
      {
        throw std::runtime_error("invalid weight in declare-weight: '" +
                                 weightToken + "'");
      }
      headPos = skipWhitespaceAndComments(smtFormula, headPos, commandEnd - 1);
      if (headPos < commandEnd - 1 && smtFormula[headPos] != ')')
      {
        throw std::runtime_error("declare-weight expects exactly one literal and weight");
      }

      result.weightDeclarations.push_back(std::move(decl));
      std::fill(result.smtWithoutProjDecls.begin() + commandStart,
                result.smtWithoutProjDecls.begin() + commandEnd,
                ' ');
    }
    pos = commandEnd;
  }
  return result;
}

std::string termSymbol(const cvc5::Term& term)
{
  return term.hasSymbol() ? term.getSymbol() : term.toString();
}

}  // namespace

TTCParser::TTCParser()
    : d_tm(),
      d_solver(ttc::createSolverWithStorage<cvc5::Solver>(d_tm)),
      d_parser(nullptr),
      d_formula(),
      d_projVars(),
      d_boolVars(),
      d_intVars(),
      d_realVars(),
      d_numConstraints(0),
      d_polytopes()
{
    d_solver.setOption("print-success", "false");
    d_solver.setOption("incremental", "true");
    d_solver.setOption("produce-models", "true");
    d_solver.setOption("produce-learned-literals", "true");
}

TTCParser::~TTCParser() {
    if (d_parser != nullptr) {
        delete d_parser;
        d_parser = nullptr;
    }
}

void TTCParser::parseFormula(const std::string &smtFormula) {
    Log(3) << "Starting parse" << std::endl;
    ProjectionDeclarationParseResult projectionDecls =
        extractExtensionDeclarations(smtFormula);
    d_hasExplicitProjectionVars = projectionDecls.sawProjectionDeclaration;
    d_hasWeights = projectionDecls.sawWeightDeclaration;
    std::istringstream iss(projectionDecls.smtWithoutProjDecls);
    d_parser = new cvc5::parser::InputParser(&d_solver);
    d_parser->setStreamInput(cvc5::modes::InputLanguage::SMT_LIB_2_6, iss, "input_stream");

    cvc5::parser::Command cmd;
    while (true) {
        cmd = d_parser->nextCommand();
        if (cmd.isNull()) {
            break;
        }
        std::stringstream out;
        cmd.invoke(&d_solver, d_parser->getSymbolManager(), out);
    }
    Log(3) << "Commands processed" << std::endl;
    auto assertions = d_solver.getAssertions();
    d_assertions = assertions;
    d_declaredVars.clear();
    if (d_parser != nullptr)
    {
        cvc5::parser::SymbolManager* sm = d_parser->getSymbolManager();
        if (sm != nullptr)
        {
            d_declaredVars = sm->getDeclaredTerms();
        }
    }
    d_projVars.clear();
    d_nonProjVars.clear();
    d_boolVars.clear();
    d_intVars.clear();
    d_realVars.clear();
    d_bvVars.clear();
    d_polytopes.clear();
    d_literalWeights.clear();
    if (d_hasExplicitProjectionVars)
    {
        std::unordered_map<std::string, cvc5::Term> declaredByName;
        for (const cvc5::Term& term : d_declaredVars)
        {
            if (term.hasSymbol())
            {
                declaredByName.emplace(term.getSymbol(), term);
            }
            declaredByName.emplace(term.toString(), term);
        }

        std::unordered_set<cvc5::Term> seenExplicitTerms;
        for (const std::string& name : projectionDecls.projectionSymbols)
        {
            auto it = declaredByName.find(name);
            if (it == declaredByName.end())
            {
                throw std::runtime_error(
                    "declare-projvar references undeclared variable '" + name + "'");
            }
            if (seenExplicitTerms.insert(it->second).second)
            {
                d_projVars.push_back(it->second);
            }
        }
    }
    if (d_hasWeights)
    {
        std::unordered_map<std::string, cvc5::Term> declaredByName;
        for (const cvc5::Term& term : d_declaredVars)
        {
            if (term.hasSymbol())
            {
                declaredByName.emplace(term.getSymbol(), term);
            }
            declaredByName.emplace(term.toString(), term);
        }

        for (const auto& decl : projectionDecls.weightDeclarations)
        {
            auto it = declaredByName.find(decl.symbol);
            if (it == declaredByName.end())
            {
                throw std::runtime_error(
                    "declare-weight references undeclared variable '" +
                    decl.symbol + "'");
            }
            if (!it->second.getSort().isBoolean())
            {
                throw std::runtime_error(
                    "declare-weight supports Boolean variables only: '" +
                    decl.symbol + "'");
            }
            LiteralWeight& weight = d_literalWeights[it->second];
            if (decl.negative)
            {
                weight.hasNegative = true;
                weight.negative = decl.weight;
                if (!weight.hasPositive)
                {
                    weight.positive = 1.0 - decl.weight;
                }
            }
            else
            {
                weight.hasPositive = true;
                weight.positive = decl.weight;
                if (!weight.hasNegative)
                {
                    weight.negative = 1.0 - decl.weight;
                }
            }
        }
    }
    if (!assertions.empty())
    {
        auto& tm = ttc::getTermBuilder(d_solver);
        cvc5::Term rawFormula = assertions.size() == 1
                                    ? assertions[0]
                                    : tm.mkTerm(cvc5::Kind::AND, assertions);
        d_formula = d_solver.simplify(rawFormula);

        // Count top-level constraints recursively.
        std::function<std::size_t(const cvc5::Term&)> countConstraints = [&](const cvc5::Term& t) {
            if (t.getKind() == cvc5::Kind::AND)
            {
                std::size_t sum = 0;
                for (size_t i = 0, n = t.getNumChildren(); i < n; ++i)
                {
                    sum += countConstraints(t[i]);
                }
                return sum;
            }
            return static_cast<std::size_t>(1);
        };
        d_numConstraints = countConstraints(d_formula);

        std::unordered_set<std::string> seen;
        std::function<void(const cvc5::Term&)> collect = [&](const cvc5::Term& t) {
            if (t.getNumChildren() == 0 && t.hasSymbol())
            {
                std::string name = termSymbol(t);
                if (name != "true" && name != "false" && seen.insert(name).second)
                {
                    if (t.getSort().isBoolean())
                    {
                        d_boolVars.insert(name);
                        Log(4) << "Boolean variable: " << name << std::endl;
                        if (!d_hasExplicitProjectionVars &&
                            name.rfind("proj_", 0) == 0)
                        {
                            d_projVars.push_back(t);
                        }
                        else
                        {
                            d_nonProjVars.push_back(t);
                        }
                    }
                    else if (t.getSort().isInteger())
                    {
                        d_intVars.insert(name);
                        d_nonProjVars.push_back(t);
                    }
                    else if (t.getSort().isReal())
                    {
                        d_realVars.insert(name);
                        d_nonProjVars.push_back(t);
                    }
                    else if (t.getSort().isBitVector())
                    {
                        d_bvVars.insert(name);
                        d_nonProjVars.push_back(t);
                    }
                }
            }
            for (size_t i = 0, n = t.getNumChildren(); i < n; ++i)
            {
                collect(t[i]);
            }
        };
        // Collect projection-variable candidates from the unsimplified formula
        // so that variables which survive in the asserted constraints but are
        // rewritten away by simplification (e.g. when the formula reduces to a
        // tautology over them) are still recognised as part of the support.
        collect(rawFormula);

        if (d_hasExplicitProjectionVars)
        {
            std::vector<cvc5::Term> filteredNonProjVars;
            filteredNonProjVars.reserve(d_nonProjVars.size());
            std::unordered_set<cvc5::Term> projectionSet(d_projVars.begin(),
                                                         d_projVars.end());
            for (const cvc5::Term& t : d_nonProjVars)
            {
                if (projectionSet.find(t) == projectionSet.end())
                {
                    filteredNonProjVars.push_back(t);
                }
            }
            d_nonProjVars.swap(filteredNonProjVars);
        }
    }
}

void TTCParser::promoteBooleanAndBvToProjection()
{
    if (d_nonProjVars.empty())
    {
        return;
    }
    std::vector<cvc5::Term> remaining;
    remaining.reserve(d_nonProjVars.size());
    for (const cvc5::Term& t : d_nonProjVars)
    {
        cvc5::Sort sort = t.getSort();
        if (sort.isBoolean() || sort.isBitVector())
        {
            d_projVars.push_back(t);
        }
        else
        {
            remaining.push_back(t);
        }
    }
    d_nonProjVars.swap(remaining);
}

#ifdef TTC_ENABLE_DDNNF
std::uint64_t TTCParser::projectedModelCount(CacheMode cacheMode,
                                             bool simplifyProp,
                                             bool contract,
                                             bool netrel,
                                             bool useAssumptions,
                                             bool useBitset,
                                             int propAt,
                                             bool useMono,
                                             bool monoTrue)
{
    Log(3) << "Starting model count with DDNNF" << std::endl;
    ProjDDNNF ddnnf(d_solver,
                    d_formula,
                    d_projVars,
                    cacheMode,
                    simplifyProp,
                    contract,
                    netrel,
                    useAssumptions,
                    useBitset,
                    propAt,
                    useMono,
                    monoTrue);
    return ddnnf.count();
}

void TTCParser::printTreeDecomposition(bool contract, bool netrel)
{
    computeProjVarOrder(d_formula, d_projVars, d_solver, true, contract, netrel);
}
#else
std::uint64_t TTCParser::projectedModelCount(CacheMode,
                                             bool,
                                             bool,
                                             bool,
                                             bool,
                                             bool,
                                             int,
                                             bool,
                                             bool)
{
    throw std::runtime_error("DDNNF-based counting is not enabled in this build");
}

void TTCParser::printTreeDecomposition(bool, bool)
{
    Log(2) << "Tree decomposition requested but DDNNF support is disabled" << std::endl;
}
#endif

// Projection-set minimization is independent of the d-DNNF/D4 stack: it only
// needs the incremental cvc5 solver, so it is always available (--arjun).
void TTCParser::minimizeProjectionSet()
{
    d_projVarsBefore = d_projVars.size();
    d_arjunWeightMultiplier = 1.0L;
    if (d_projVars.empty())
    {
        d_projVarsAfter = 0;
        return;
    }
    Log(3) << "Minimizing projection set of size " << d_projVarsBefore
           << std::endl;

    arjun::Reduction r =
        arjun::analyzeProjectionSet(d_solver, d_formula, d_projVars);

    auto weightOf = [&](const cvc5::Term& v, bool value) -> long double {
        auto it = d_literalWeights.find(v);
        LiteralWeight w =
            (it != d_literalWeights.end()) ? it->second : LiteralWeight{};
        return static_cast<long double>(value ? w.positive : w.negative);
    };

    std::vector<cvc5::Term> kept = std::move(r.support);

    // Backbone variables are constant in every model. Dropping them never
    // changes the unweighted projected count; under weighted counting each
    // contributes the fixed weight of its forced value.
    for (const auto& [v, value] : r.forced)
    {
        if (d_hasWeights)
        {
            d_arjunWeightMultiplier *= weightOf(v, value);
        }
        Log(4) << "  forced " << v << " = " << (value ? "true" : "false")
               << std::endl;
    }

    // Implicitly-defined variables are functionally determined by the support,
    // so dropping one preserves the unweighted count. Under weighted counting
    // the removed variable's value -- and thus its weight factor -- varies per
    // model, so removal is only sound when its two weights are equal (the
    // factor is then constant); otherwise it must stay in the projection set.
    std::size_t definedDropped = 0;
    for (const cvc5::Term& v : r.defined)
    {
        if (!d_hasWeights)
        {
            ++definedDropped;
            continue;
        }
        auto it = d_literalWeights.find(v);
        LiteralWeight w =
            (it != d_literalWeights.end()) ? it->second : LiteralWeight{};
        if (w.positive == w.negative)
        {
            d_arjunWeightMultiplier *= static_cast<long double>(w.positive);
            ++definedDropped;
        }
        else
        {
            kept.push_back(v);  // unsafe to drop under asymmetric weights
        }
    }

    // Unconstrained variables occur in no constraint, so each one doubles the
    // model count: dropping it from the projection set requires multiplying the
    // unweighted count back by 2, and the weighted count by the sum of its two
    // literal weights (the contribution of summing over both values).
    for (const cvc5::Term& v : r.free)
    {
        if (d_hasWeights)
        {
            auto it = d_literalWeights.find(v);
            LiteralWeight w =
                (it != d_literalWeights.end()) ? it->second : LiteralWeight{};
            d_arjunWeightMultiplier *=
                static_cast<long double>(w.positive + w.negative);
        }
        else
        {
            // The model count gains one factor per value the variable can take:
            // 2 for a Boolean, 2^width for a bit-vector.
            cvc5::Sort sort = v.getSort();
            unsigned width = sort.isBitVector() ? sort.getBitVectorSize() : 1u;
            // A 64-bit unsigned count cannot represent 2^64 or beyond; leave such
            // (degenerate) variables in the projection set rather than overflow.
            if (width < 64u)
            {
                d_arjunCountMultiplier *= (1ULL << width);
            }
            else
            {
                kept.push_back(v);
            }
        }
        Log(4) << "  unconstrained " << v << std::endl;
    }

    d_projVars = std::move(kept);
    d_projVarsAfter = d_projVars.size();
    Log(3) << "Projection set reduced to " << d_projVarsAfter << " (forced "
           << r.forced.size() << ", defined " << definedDropped
           << ", unconstrained " << r.free.size() << ")" << std::endl;
}

std::optional<bool> TTCParser::checkSatWithSTP(const std::string& smtFormula)
{
    // Require only Boolean and bit-vector variables.
    if (!d_intVars.empty() || !d_realVars.empty())
    {
        return std::nullopt;
    }

    return csb::stpCheckSat(smtFormula);
}

void TTCParser::computePolytopes()
{
    d_polytopes.clear();
    if (d_formula.isNull())
    {
        return;
    }

    try
    {
        auto abstraction = ttc::getBooleanAbstractionAig(d_solver, d_formula);
        auto result = ttc::enumeratePolytopes(abstraction, d_solver);
        d_polytopes = std::move(result.polytopes);
    }
    catch (const std::exception&)
    {
        try
        {
            auto manual = ttc::buildBooleanAigFromTerm(d_formula);
            auto result = ttc::enumeratePolytopes(manual, d_solver);
            d_polytopes = std::move(result.polytopes);
        }
        catch (const std::exception&)
        {
            d_polytopes.clear();
        }
    }
}

std::vector<cvc5::Term> TTCParser::realVariables() const
{
    std::unordered_set<cvc5::Term> seen;
    std::vector<cvc5::Term> vars;
    vars.reserve(d_nonProjVars.size() + d_projVars.size());
    for (const auto& term : d_nonProjVars)
    {
        if (term.getSort().isReal() && seen.insert(term).second)
        {
            vars.push_back(term);
        }
    }
    for (const auto& term : d_projVars)
    {
        if (term.getSort().isReal() && seen.insert(term).second)
        {
            vars.push_back(term);
        }
    }
    return vars;
}
