#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <optional>
#include <utility>
#include <cctype>
#include <unistd.h>

#include <boost/multiprecision/cpp_int.hpp>

#include "cache_mode.hpp"
#if defined(TTC_ENABLE_DDNNF)
#include "d4_interface.hpp"
#endif
#include "features.hpp"
#include "logger.hpp"
#include "pact/pact.hpp"
#include "projected_enumerator.hpp"
#include "parser.hpp"
#include "profiler.hpp"
#include "tocnf/tocnf.hpp"
#include "eager/bvcnf.hpp"
#include "weighted/tounweighted.hpp"
#if defined(TTC_ENABLE_DDNNF)
#include "var_order.hpp"
#endif
#include "volume/lra_volume.hpp"

namespace po = boost::program_options;

namespace
{
using boost::multiprecision::cpp_int;

bool g_useNativeXor = false;
bool g_hasCadicalXorSupport = false;
bool g_xorTraceEnabled = false;

void printWeightedCount(long double value)
{
  std::ios::fmtflags oldFlags = std::cout.flags();
  std::streamsize oldPrecision = std::cout.precision();
  std::cout.unsetf(std::ios::floatfield);
  std::cout << std::setprecision(15) << static_cast<double>(value);
  std::cout.flags(oldFlags);
  std::cout.precision(oldPrecision);
}

void markPactRunningInBanner()
{
  std::string info = banner_solver_info();
  const std::string marker = "counter: pact is running";
  if (info.find(marker) == std::string::npos)
  {
    if (!info.empty())
    {
      info += "; ";
    }
    info += marker;
    set_banner_solver_info(std::move(info));
  }
}

template <typename SolverT>
auto detectCadicalXorSupportImpl(SolverT& solver, int)
    -> decltype(solver.hasCadicalXorSupport(), bool())
{
  return solver.hasCadicalXorSupport();
}

template <typename SolverT>
bool detectCadicalXorSupportImpl(SolverT& solver, long)
{
  try
  {
    solver.setOption("sat-use-native-xor", "true");
    return true;
  }
  catch (const cvc5::CVC5ApiException&)
  {
    return false;
  }
}

template <typename SolverT>
bool detectCadicalXorSupport(SolverT& solver)
{
  bool hasSupport = detectCadicalXorSupportImpl(solver, 0);
  if (!hasSupport)
  {
    hasSupport = detectCadicalXorSupportImpl(solver, 0L);
  }
  return hasSupport;
}

template <typename SolverT>
auto applyNativeXorImpl(SolverT& solver, bool enable, int)
    -> decltype(solver.setSatUseNativeXor(enable), void())
{
  solver.setSatUseNativeXor(enable);
}

template <typename SolverT>
void applyNativeXorImpl(SolverT& solver, bool enable, long)
{
  solver.setOption("sat-use-native-xor", enable ? "true" : "false");
}

template <typename SolverT>
void applyNativeXor(SolverT& solver, bool enable)
{
  applyNativeXorImpl(solver, enable, 0);
}

std::optional<cpp_int> computeTotalAssignments(
    const std::vector<cvc5::Term>& projectionVars)
{
  cpp_int total = 1;
  for (const cvc5::Term& var : projectionVars)
  {
    cvc5::Sort sort = var.getSort();
    if (sort.isBoolean())
    {
      total <<= 1;
    }
    else if (sort.isBitVector())
    {
      total <<= sort.getBitVectorSize();
    }
    else
    {
      return std::nullopt;
    }
  }
  return total;
}

struct WeightedPbResult
{
  long double weightedCount = 0.0L;
  std::uint64_t unweightedCount = 0;
  std::uint64_t smtCalls = 0;
  int originalVars = 0;
  int originalClauses = 0;
  int convertedVars = 0;
  int convertedClauses = 0;
  int samplingVars = 0;
  int divisorPower = 0;
  long double multiplier = 1.0L;
};

int findOrAddCnfVar(ToCNF::CNFFormula& cnf, const cvc5::Term& term)
{
  auto termIt = cnf.termToIdx.find(term);
  if (termIt != cnf.termToIdx.end())
  {
    return termIt->second;
  }

  auto strIt = cnf.varToIdx.find(term.toString());
  if (strIt != cnf.varToIdx.end())
  {
    cnf.termToIdx.emplace(term, strIt->second);
    return strIt->second;
  }

  if (term.hasSymbol())
  {
    auto symIt = cnf.varToIdx.find(term.getSymbol());
    if (symIt != cnf.varToIdx.end())
    {
      cnf.termToIdx.emplace(term, symIt->second);
      return symIt->second;
    }
  }

  int id = ++cnf.varCount;
  if (static_cast<int>(cnf.idxToTerm.size()) <= id)
  {
    cnf.idxToTerm.resize(id + 1);
  }
  cnf.idxToTerm[id] = term;
  cnf.termToIdx.emplace(term, id);
  cnf.varToIdx.emplace(term.toString(), id);
  if (term.hasSymbol())
  {
    cnf.varToIdx.emplace(term.getSymbol(), id);
  }
  return id;
}

void configureWeightedPbSolver(cvc5::Solver& solver,
                               bool useBvPact,
                               bool useNativeXor)
{
  try { solver.setOption("print-success", "false"); }
  catch (const cvc5::CVC5ApiException&) {}
  try { solver.setOption("incremental", "true"); }
  catch (const cvc5::CVC5ApiException&) {}
  try { solver.setOption("produce-models", "true"); }
  catch (const cvc5::CVC5ApiException&) {}
  try { solver.setOption("produce-learned-literals", "true"); }
  catch (const cvc5::CVC5ApiException&) {}
  try { solver.setOption("bv-sat-solver", "cryptominisat"); }
  catch (const cvc5::CVC5ApiException&) {}
  if (useBvPact)
  {
    try { solver.setOption("bv-to-bool", "false"); }
    catch (const cvc5::CVC5ApiException&) {}
  }
  try { solver.setLogic("QF_BV"); }
  catch (const cvc5::CVC5ApiException&) {}
  try { applyNativeXor(solver, useNativeXor); }
  catch (const cvc5::CVC5ApiException&) {}
}

void assertCnf(cvc5::Solver& solver,
               const std::vector<cvc5::Term>& variables,
               const std::vector<std::vector<int>>& clauses)
{
  auto& tm = ttc::getTermBuilder(solver);
  for (const auto& clause : clauses)
  {
    if (clause.empty())
    {
      solver.assertFormula(tm.mkBoolean(false));
      continue;
    }

    std::vector<cvc5::Term> terms;
    terms.reserve(clause.size());
    for (int lit : clause)
    {
      int var = std::abs(lit);
      if (var <= 0 || var >= static_cast<int>(variables.size()))
      {
        throw std::runtime_error("CNF literal is outside the variable range");
      }
      cvc5::Term term = variables[var];
      if (lit < 0)
      {
        term = tm.mkTerm(cvc5::Kind::NOT, {term});
      }
      terms.push_back(term);
    }

    if (terms.size() == 1)
    {
      solver.assertFormula(terms[0]);
    }
    else
    {
      solver.assertFormula(tm.mkTerm(cvc5::Kind::OR, terms));
    }
  }
}

WeightedPbResult runWeightedProjectionBoost(TTCParser& parser,
                                            std::uint64_t seed,
                                            bool useNativeXor,
                                            bool useBvPact,
                                            double epsilon = 0.8,
                                            double delta = 0.2)
{
  for (const cvc5::Term& term : parser.projectionVars())
  {
    if (!term.getSort().isBoolean())
    {
      throw std::runtime_error(
          "weighted --PB currently supports Boolean projection variables only");
    }
  }

  ToCNF cnfBuilder(parser.solver(), parser.assertions());
  ToCNF::CNFFormula cnf = cnfBuilder.build();

  ttc::weighted::WeightedCnf weighted;
  weighted.clauses = cnf.clauses;
  std::unordered_set<int> seenSampling;
  const auto& weights = parser.literalWeights();
  for (const cvc5::Term& term : parser.projectionVars())
  {
    int id = findOrAddCnfVar(cnf, term);
    if (seenSampling.insert(id).second)
    {
      weighted.samplingVars.push_back(id);
    }

    TTCParser::LiteralWeight literalWeight;
    auto weightIt = weights.find(term);
    if (weightIt != weights.end())
    {
      literalWeight = weightIt->second;
    }
    weighted.weights.push_back(
        {id,
         static_cast<long double>(literalWeight.positive),
         static_cast<long double>(literalWeight.negative)});
  }
  weighted.varCount = cnf.varCount;

  ttc::weighted::ToUnweightedConverter converter;
  ttc::weighted::UnweightedCnf unweighted = converter.convert(weighted);

  ttc::TermBuilderHelper<cvc5::Solver>::storage_type storage;
  cvc5::Solver solver = ttc::createSolverWithStorage<cvc5::Solver>(storage);
  configureWeightedPbSolver(solver, useBvPact, useNativeXor);
  auto& tm = ttc::getTermBuilder(solver);
  cvc5::Sort boolSort = tm.getBooleanSort();
  std::vector<cvc5::Term> cnfVars(unweighted.varCount + 1);
  for (int i = 1; i <= unweighted.varCount; ++i)
  {
    cnfVars[i] = tm.mkConst(boolSort, "__ttc_w2u_" + std::to_string(i));
  }
  assertCnf(solver, cnfVars, unweighted.clauses);

  std::vector<cvc5::Term> projectionVars;
  projectionVars.reserve(unweighted.samplingVars.size());
  for (int var : unweighted.samplingVars)
  {
    if (var <= 0 || var > unweighted.varCount)
    {
      throw std::runtime_error("converted sampling variable is out of range");
    }
    projectionVars.push_back(cnfVars[var]);
  }

  Pact counter(
      solver, projectionVars, seed, useNativeXor, useBvPact, epsilon, delta);
  std::uint64_t unweightedCount = counter.count();

  WeightedPbResult result;
  result.unweightedCount = unweightedCount;
  result.smtCalls = counter.getSmtCallCount();
  result.originalVars = unweighted.originalVarCount;
  result.originalClauses = unweighted.originalClauseCount;
  result.convertedVars = unweighted.varCount;
  result.convertedClauses = static_cast<int>(unweighted.clauses.size());
  result.samplingVars = static_cast<int>(unweighted.samplingVars.size());
  result.divisorPower = unweighted.divisorPower;
  result.multiplier = unweighted.multiplier;
  result.weightedCount =
      std::ldexp(static_cast<long double>(unweightedCount) *
                     unweighted.multiplier,
                 -unweighted.divisorPower);
  return result;
}

cpp_int clampNonNegative(cpp_int value)
{
  if (value < 0)
  {
    return cpp_int(0);
  }
  return value;
}

std::optional<std::uint64_t> countUnsatAssignmentsQuantified(
    TTCParser& parser, std::uint64_t& smtCallsOut, std::uint64_t seed)
{
  cvc5::Solver& solver = parser.solver();
  const std::vector<cvc5::Term>& nonProj = parser.nonProjectionVars();
  const std::vector<cvc5::Term>& proj = parser.projectionVars();
  const cvc5::Term& baseFormula = parser.formula();
  if (baseFormula.isNull())
  {
    return std::nullopt;
  }

  auto& tm = ttc::getTermBuilder(solver);
  cvc5::Term unsatTerm;
  if (nonProj.empty())
  {
    unsatTerm = tm.mkTerm(cvc5::Kind::NOT, {baseFormula});
  }
  else
  {
    std::vector<cvc5::Term> boundVars;
    boundVars.reserve(nonProj.size());
    for (std::size_t i = 0, n = nonProj.size(); i < n; ++i)
    {
      std::ostringstream name;
      name << "unsat_q" << i;
      boundVars.push_back(tm.mkVar(nonProj[i].getSort(), name.str()));
    }
    cvc5::Term substituted =
        baseFormula.substitute(nonProj, boundVars);
    cvc5::Term boundList = tm.mkTerm(cvc5::Kind::VARIABLE_LIST, boundVars);
    cvc5::Term existsTerm =
        tm.mkTerm(cvc5::Kind::EXISTS, {boundList, substituted});
    cvc5::Term quantified = tm.mkTerm(cvc5::Kind::NOT, {existsTerm});
    try
    {
      unsatTerm = solver.getQuantifierElimination(quantified);
    }
    catch (const cvc5::CVC5ApiException&)
    {
      unsatTerm = quantified;
    }
  }

  cvc5::Solver quantSolver = ttc::makeSolverWithBuilder(solver);
  quantSolver.setOption("print-success", "false");
  quantSolver.setOption("incremental", "true");
  quantSolver.setOption("produce-models", "true");
  if (g_hasCadicalXorSupport)
  {
    applyNativeXor(quantSolver, g_useNativeXor);
  }
  // quantSolver.setXorAssertionVerbose(g_xorTraceEnabled);
  if (!quantSolver.isLogicSet())
  {
    quantSolver.setLogic("ALL");
  }
  quantSolver.assertFormula(unsatTerm);
  Pact unsatCounter(quantSolver, proj, seed, g_useNativeXor);
  std::uint64_t unsatCount = unsatCounter.count();
  smtCallsOut += unsatCounter.getSmtCallCount();
  return unsatCount;
}

std::optional<std::uint64_t> countUnsatAssignmentsQuantified(
    TTCParser& parser, std::uint64_t seed)
{
  std::uint64_t dummy = 0;
  return countUnsatAssignmentsQuantified(parser, dummy, seed);
}

std::string overrideLogicForQuantifiers(std::string input)
{
  const std::string directive = "(set-logic";
  bool replaced = false;
  std::size_t pos = 0;
  while ((pos = input.find(directive, pos)) != std::string::npos)
  {
    std::size_t start = pos + directive.size();
    start = input.find_first_not_of(" \t\r\n", start);
    if (start == std::string::npos)
    {
      break;
    }
    std::size_t end = input.find_first_of(" \t\r\n)", start);
    if (end == std::string::npos)
    {
      break;
    }
    input.replace(start, end - start, "ALL");
    pos = start + 3;
    replaced = true;
  }
  if (!replaced)
  {
    input.insert(0, "(set-logic ALL)\n");
  }
  return input;
}

// Decide the default --PB backend by peeking at the input file. cvcxor
// (Gauss-Jordan in cvc5's propagator) is much faster on LRA/theory inputs, but
// it is NOT a safe blanket default:
//   * bit-vector inputs are far better served by BV-PACT (bit-blast XOR into
//     the BV SAT solver); cvcxor stalls on them (e.g. gaussj).
//   * cvcxor runs the count solver with simplification=none, which does not
//     inline define-fun bodies, so inputs that use define-fun error out
//     ("Function terms are only supported with higher-order logic").
// We therefore default to cvcxor only when the input has no bit-vectors and no
// define-fun, and fall back to BV-PACT otherwise. --cvcxor / --bv-pact still
// force either backend explicitly.
bool preferCvcXorDefault(const std::string& filename)
{
  std::ifstream in(filename);
  if (!in)
  {
    return false;
  }
  std::string input((std::istreambuf_iterator<char>(in)),
                    std::istreambuf_iterator<char>());
  // define-fun (function terms) are incompatible with simplification=none.
  if (input.find("define-fun") != std::string::npos)
  {
    return false;
  }
  // Inspect the (set-logic ...) directive for bit-vectors.
  const std::string directive = "(set-logic";
  std::size_t pos = input.find(directive);
  if (pos != std::string::npos)
  {
    std::size_t start = pos + directive.size();
    start = input.find_first_not_of(" \t\r\n", start);
    if (start != std::string::npos)
    {
      std::size_t end = input.find_first_of(" \t\r\n)", start);
      if (end != std::string::npos)
      {
        std::string logic = input.substr(start, end - start);
        std::transform(logic.begin(), logic.end(), logic.begin(),
                       [](unsigned char ch) {
                         return static_cast<char>(std::toupper(ch));
                       });
        // "ALL" may admit bit-vectors: be conservative and prefer BV-PACT.
        if (logic == "ALL" || logic.find("BV") != std::string::npos)
        {
          return false;
        }
      }
    }
  }
  return true;
}

std::string ensureLogicIncludesBv(std::string input)
{
  const std::string directive = "(set-logic";
  std::size_t pos = 0;
  while ((pos = input.find(directive, pos)) != std::string::npos)
  {
    std::size_t start = pos + directive.size();
    start = input.find_first_not_of(" \t\r\n", start);
    if (start == std::string::npos)
    {
      break;
    }
    std::size_t end = input.find_first_of(" \t\r\n)", start);
    if (end == std::string::npos)
    {
      break;
    }
    std::string logic = input.substr(start, end - start);
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
      input.replace(start, end - start, newLogic);
      pos = start + newLogic.size();
    }
    else
    {
      pos = end;
    }
  }
  return input;
}

// ---------------------------------------------------------------------------
// Automatic engine dispatch helpers.
//
// ttc selects one of three counting engines from the (set-logic ...) directive
// and the projection variables of the input:
//   1. logic BV / UFBV                          -> eager bit-blast + ApproxMC
//   2. other theories, BV/Bool projection vars  -> projection counting (--PB)
//   3. QF_LRA with no projection variables       -> LRA volume computation
//   4. otherwise                                 -> unsupported (reported)
// ---------------------------------------------------------------------------

// Read the upper-cased logic token from the input's (set-logic ...) directive.
// Returns the empty string when no directive is present.
std::string extractLogicUpper(const std::string& filename)
{
  std::ifstream in(filename);
  if (!in)
  {
    return std::string();
  }
  std::string input((std::istreambuf_iterator<char>(in)),
                    std::istreambuf_iterator<char>());
  const std::string directive = "(set-logic";
  std::size_t pos = input.find(directive);
  if (pos == std::string::npos)
  {
    return std::string();
  }
  std::size_t start = input.find_first_not_of(" \t\r\n", pos + directive.size());
  if (start == std::string::npos)
  {
    return std::string();
  }
  std::size_t end = input.find_first_of(" \t\r\n)", start);
  if (end == std::string::npos)
  {
    return std::string();
  }
  std::string logic = input.substr(start, end - start);
  std::transform(logic.begin(), logic.end(), logic.begin(),
                 [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
  return logic;
}

// Remove every occurrence of token from s (in place).
void eraseLogicToken(std::string& s, const std::string& token)
{
  std::size_t pos;
  while ((pos = s.find(token)) != std::string::npos)
  {
    s.erase(pos, token.size());
  }
}

// True when the logic is a pure bit-vector logic (BV / UFBV / QF_BV / QF_UFBV);
// arrays, floating point, arithmetic, etc. are NOT pure-BV.
bool logicIsPureBv(const std::string& logicUpper)
{
  if (logicUpper.empty())
  {
    return false;
  }
  std::string t = logicUpper;
  if (t.rfind("QF_", 0) == 0)
  {
    t = t.substr(3);
  }
  eraseLogicToken(t, "UF");
  bool hadBv = t.find("BV") != std::string::npos;
  eraseLogicToken(t, "BV");
  return hadBv && t.empty();
}

// True when every projection variable is Boolean or a bit-vector.
bool projVarsAreBvOrBool(const std::vector<cvc5::Term>& projVars)
{
  for (const cvc5::Term& v : projVars)
  {
    cvc5::Sort sort = v.getSort();
    if (!sort.isBoolean() && !sort.isBitVector())
    {
      return false;
    }
  }
  return true;
}

// Build a writable path for the intermediate model-preserving CNF that the
// bit-vector ApproxMC path produces under automatic dispatch (the user did not
// supply one via --bvcnf <file>).
std::string makeAutoCnfPath(const std::string& filename)
{
  const char* tmp = std::getenv("TMPDIR");
  std::string dir = (tmp != nullptr && *tmp != '\0') ? tmp : "/tmp";
  std::string base = filename;
  std::size_t slash = base.find_last_of("/\\");
  if (slash != std::string::npos)
  {
    base = base.substr(slash + 1);
  }
  if (base.empty())
  {
    base = "ttc";
  }
  return dir + "/" + base + "." + std::to_string(static_cast<long>(::getpid())) +
         ".bvcnf";
}

// Report that the input matches none of the auto-dispatch engine classes and
// explain how the engine is chosen.
void printUnsupportedEngine(const std::string& logicUpper,
                            std::size_t numProjVars)
{
  std::cerr
      << "Error: unsupported input for automatic engine selection.\n"
      << "  logic                : "
      << (logicUpper.empty() ? "(none)" : logicUpper) << "\n"
      << "  projection variables : " << numProjVars << "\n\n"
      << "ttc selects a counting engine from the logic and projection variables:\n"
      << "  1. logic BV / UFBV                          -> bit-vector ApproxMC "
         "(eager bit-blast)\n"
      << "  2. other theories with BV/Bool projection   -> projection counting "
         "(--pact)\n"
      << "  3. QF_LRA with no projection variables       -> LRA volume "
         "computation\n\n"
      << "This input matches none of these. Declare BV/Bool projection "
         "variables,\n"
      << "use a supported logic, or force an engine explicitly (e.g. --pact, "
         "--bvcnf, --enum)."
      << std::endl;
}

}  // namespace

int main(int argc, char *argv[]) {
  // General options shared by every engine. By default ttc auto-selects an
  // engine from the (set-logic ...) directive and the projection variables;
  // the options below are grouped by the engine they configure.
  po::options_description general("General options");
  general.add_options()
      ("help,h", "Display help information")
      ("approx,a",
       "Use approximate model counting")
      ("seed,s",
       po::value<std::uint64_t>()->default_value(42),
       "Seed for randomized approximate counting")
      ("epsilon,e",
       po::value<double>()->default_value(0.8, "0.8"),
       "Approximation tolerance: the count is within a (1+epsilon) "
       "multiplicative factor of the true count (ApproxMC epsilon)")
      ("delta,d",
       po::value<double>()->default_value(0.2, "0.2"),
       "Confidence: the (1+epsilon) guarantee holds with probability at least "
       "1-delta (ApproxMC delta)")
      ("verbose,v", po::value<int>()->default_value(1),
       "Verbosity level (0-10)")
      ("verbosity-time,vt",
       po::value<double>()->default_value(5.0),
       "Print progress every K seconds")
      ("trace,t",
       po::value<std::vector<std::string>>()->composing()->multitoken(),
       "Enable tracing for selected components");

  // Engine 1: pure bit-vector logics (BV / UFBV) -> eager bit-blast + ApproxMC.
  po::options_description classBv(
      "Engine 1 -- bit-vector counting (auto: logic BV / UFBV)");
  classBv.add_options()
      ("bvcnf", po::value<std::string>(),
       "Eagerly bit-blast a QF_BV formula to a model-preserving CNF, write it "
       "to the given file, and run ApproxMC to count its bit-vector models "
       "(forces this engine; auto-selected for BV/UFBV logics)");

  // Engine 2: other theories with BV/Bool projection variables -> --pact.
  po::options_description classPb(
      "Engine 2 -- projection counting (auto: other theories, BV/Bool "
      "projection vars)");
  classPb.add_options()
      ("pact,P",
       "Enable projection-based approximate counting (forces this engine; "
       "auto-selected for non-BV/non-LRA logics with BV/Bool projection "
       "variables)")
      ("enum,E", "Enumerate projected assignments exactly")
      ("xor", po::value<std::string>(),
       "XOR (parity hash) backend: 'cvc' = Gauss-Jordan inside cvc5's "
       "propagator (default mode); 'cadical' = CaDiCaL native XOR engine; "
       "'cms' = route to the bit-vector SAT solver (CryptoMiniSat); "
       "'blast' = Tseitin-blast to CNF for plain CaDiCaL. Unset: 'cvc' for "
       "LRA/theory inputs, auto-fallback to 'cms' for bit-vector inputs.")
      ("xor-activation", po::value<std::string>(),
       "How the native-XOR galloping search switches hashes on/off: 'rebuild' "
       "(default) rebuilds the count solver per galloping level so it sees only "
       "the active hashes; 'literal' folds an indicator variable into each "
       "parity and toggles hashes by assumption on one solver. 'rebuild' is "
       "measured faster here -- 'literal' accumulates every explored hash on "
       "the solver, which outweighs the saved rebuild.");

  // Engine 3: pure QF_LRA with no projection variables -> LRA volume.
  po::options_description classLra(
      "Engine 3 -- LRA volume computation (auto: QF_LRA, no projection vars; "
      "fully automatic, no tunable options)");

  // Engine 4 (d-DNNF / exact CNF counting) and its tuning knobs are inactive by
  // default and hidden from --help. They are still accepted on the command line
  // so existing invocations keep working.
  po::options_description advanced("Advanced / exact counting options");
  advanced.add_options()
      ("arjun", "Minimize projection set before counting")
      ("unsat",
       "Count projection assignments that make the constraints unsatisfiable")
      ("unsat-q",
       "Use quantified reasoning when counting unsatisfied projection assignments")
      ("rewrite-prop", "Use rewrite-based propagation instead of cvc5 simplify")
      ("no-contract",
       "Do not contract to projection variables for tree decomposition")
      ("netrel", "Apply netlist relevance reduction before contraction")
      ("tocnf",
       "Convert assertions to CNF")
#if defined(TTC_ENABLE_DDNNF)
      ("d4,D", "Use D4 for model counting on CNF abstraction")
      ("full-decomp",
       "Do not stop D4 decomposition due to SMT dependencies")
#endif
      ("output,o", po::value<std::string>(),
       "CNF output file")
      ("FB", po::value<std::string>(),
       "Write boolean abstraction (CNF and mapping) to file")
      ("cache", po::value<std::string>()->default_value("synt"),
       "Cache hit strategy: synt, sem, hash, canon, bool")
#if defined(TTC_ENABLE_DDNNF)
      ("smt-cache", "Enable SMT-based component caching for D4")
      ("no-decompose",
       "Disable formula decomposition into independent components")
#endif
      ("assumptions", "Use incremental assumptions instead of push/pop")
      ("bitset-vars", "Use bitset-based variable tracking for decomposition")
      ("propat", po::value<int>()->default_value(1),
       "Propagate only every K assignments")
      ("mono", "Enable monotonic reasoning with false as default assignment")
      ("mono-true",
       "Enable monotonic reasoning with true as default assignment")
#if defined(TTC_ENABLE_DDNNF)
      ("no-res-simp",
       "Disable residual SMT simplification before D4 decomposition")
#endif
      ;

  po::options_description hidden;
  hidden.add_options()
      ("input-file", po::value<std::string>(), "SMT2 file");

  // Visible help groups the options by the three active engines; the advanced
  // (engine 4 / exact CNF) options and the positional input-file are accepted
  // by the parser but hidden from --help.
  po::options_description visible("Allowed options");
  visible.add(general).add(classBv).add(classPb).add(classLra);
  po::options_description all;
  all.add(visible).add(advanced).add(hidden);

  po::positional_options_description p;
  p.add("input-file", 1);

  po::variables_map vm;
  try {
    std::vector<std::string> args;
    args.reserve(argc);
    for (int i = 1; i < argc; ++i) {
      std::string a(argv[i]);
      if (a == "-vt")
        a = "--verbosity-time";
      else if (a == "-PB" || a == "--PB")  // legacy aliases for --pact
        a = "--pact";
      else if (a == "-PE" || a == "--PE")  // legacy aliases for --enum
        a = "--enum";
      args.push_back(a);
    }
    po::store(po::command_line_parser(args).options(all).positional(p).run(),
              vm);
    po::notify(vm);
  } catch (const std::exception &ex) {
    std::cerr << "Error parsing command line: " << ex.what() << std::endl;
    return 1;
  }

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << "Usage: " << argv[0] << " [options] filename.smt2"
              << std::endl;
    std::cout << visible << std::endl;
    std::cout
        << "Automatic engine selection (when no engine is forced):\n"
        << "  1. logic BV / UFBV                          -> bit-vector "
           "ApproxMC (eager bit-blast)\n"
        << "  2. other theories with BV/Bool projection   -> projection "
           "counting (--pact)\n"
        << "  3. QF_LRA with no projection variables       -> LRA volume "
           "computation\n"
        << "  4. otherwise                                 -> unsupported "
           "(reports an error and exits)\n"
        << std::endl;
    return 1;
  }

  std::string filename = vm["input-file"].as<std::string>();
  bool simplifyProp = vm.count("rewrite-prop") == 0;
  bool noContract = vm.count("no-contract") > 0;
  bool useNetrel = vm.count("netrel") > 0;
  bool useArjun = vm.count("arjun") > 0;
  bool useApprox = vm.count("approx") > 0;
  bool useProjectionBoost = vm.count("pact") > 0;
  bool useProjectionEnumerate = vm.count("enum") > 0;
  std::uint64_t approxSeed = vm["seed"].as<std::uint64_t>();
  double approxEpsilon = vm["epsilon"].as<double>();
  double approxDelta = vm["delta"].as<double>();
  bool useUnsatQuant = vm.count("unsat-q") > 0;
  bool countUnsatAssignments = vm.count("unsat") > 0 || useUnsatQuant;

  // XOR (parity hash) backend for projection counting:
  //   --xor cvc      -> Gauss-Jordan inside cvc5's CadicalPropagator (a single
  //                     propagator does theory + parity, pruning the search
  //                     like CMS/ApproxMC). The current default mode.
  //   --xor cadical  -> hand each parity to CaDiCaL's native XOR engine.
  //   --xor cms      -> route the XOR hashes to the bit-vector SAT solver
  //                     (CryptoMiniSat); the BV-PACT substitution path.
  //   --xor blast    -> Tseitin-blast the XOR hashes to CNF and solve them as
  //                     plain clauses with CaDiCaL (no Gauss-Jordan).
  //   (unset)        -> smart default: 'cvc' for LRA/theory inputs, auto-falls
  //                     back to 'cms' for bit-vector / define-fun inputs where
  //                     Gauss-Jordan stalls or cannot run.
  std::string xorMode;
  if (vm.count("xor"))
  {
    xorMode = vm["xor"].as<std::string>();
    if (xorMode != "cvc" && xorMode != "cadical" && xorMode != "cms"
        && xorMode != "blast")
    {
      std::cerr << "Error: --xor expects 'cvc', 'cadical', 'cms' or 'blast'"
                << std::endl;
      return 1;
    }
  }
  // --xor-activation {literal,rebuild}: how the native-XOR galloping search
  // toggles hashes. Default 'literal' (indicator-assumption, no per-level
  // rebuild). Carried into Pact via TTC_XOR_ACTIVATION.
  if (vm.count("xor-activation"))
  {
    std::string act = vm["xor-activation"].as<std::string>();
    if (act != "literal" && act != "rebuild")
    {
      std::cerr << "Error: --xor-activation expects 'literal' or 'rebuild'"
                << std::endl;
      return 1;
    }
    setenv("TTC_XOR_ACTIVATION", act.c_str(), 1);
  }

  const bool xorForceCvc = (xorMode == "cvc");
  const bool xorForceCadical = (xorMode == "cadical");
  const bool xorForceCms = (xorMode == "cms");
  bool useCvcXor = false;
  bool useXorCnf = (xorMode == "blast");
  bool requestNativeXor = false;

  // ----- Automatic engine dispatch (pre-parse classification) ---------------
  // When the user does not force an engine, classify the input from its
  // (set-logic ...) directive. Projection-variable types are only known after
  // parsing, so pure-LRA is resolved into PB vs. volume later; everything that
  // can run projection counting is configured here so the --PB backend (cvcxor
  // / bv-pact) is set up exactly as if --PB had been passed explicitly.
  const bool explicitEngine =
      vm.count("pact") || vm.count("enum") || vm.count("bvcnf")
      || vm.count("tocnf") || vm.count("d4") || vm.count("FB")
      || vm.count("unsat") || vm.count("unsat-q");
  const bool autoDispatch = !explicitEngine;
  bool autoBvCnf = false;     // rule 1: BV/UFBV -> eager bit-blast + ApproxMC
  bool autoDeferred = false;  // rule 2 (PB) or rule 3 (volume), resolved after parsing
  std::string autoLogic;
  if (autoDispatch)
  {
    autoLogic = extractLogicUpper(filename);
    if (logicIsPureBv(autoLogic))
    {
      autoBvCnf = true;
    }
    else
    {
      // Rule 2 (projection counting) or rule 3 (LRA volume): the choice depends
      // on the projection-variable sorts, which are only known after parsing.
      // The PB backend (cvcxor / bv-pact) is configured below so that, if this
      // resolves to projection counting, it behaves exactly like --PB.
      autoDeferred = true;
    }
  }

  // The XOR backend is configured for the deferred case too so that, if it
  // resolves to projection counting, it behaves exactly like --pact. These
  // settings are inert for the volume path.
  const bool pbBackendNeeded = useProjectionBoost || autoDeferred;
  if (xorForceCvc)
  {
    useCvcXor = true;
  }
  else if (xorForceCadical)
  {
    // Hand each parity to CaDiCaL's native XOR engine (no substitution path).
    requestNativeXor = true;
    // The new ~/solvers/cadical ships an incremental Gauss-Jordan XOR engine
    // (a watched-column port of CMS's EGaussian) but defaults it OFF
    // (gauss=0) with CNF-blasting ON (xorblast=1). With those defaults every
    // parity hash is Tseitin-expanded to plain CNF and solved with no Gaussian
    // reasoning -- exponential on the dense XOR hashes, so counting hangs at the
    // first hash level. Enable the GJ engine exactly as ApproxMC drives this
    // same build (gauss=1, xorblast=0); CaDiCaL reads these per-option
    // CADICAL_<NAME> env vars when each solver is constructed. (The --xor cvc
    // path solves parities in cvc5's own propagator instead, so it leaves
    // CaDiCaL's engine at its defaults.)
    setenv("CADICAL_GAUSS", "1", 1);
    setenv("CADICAL_XORBLAST", "0", 1);
  }
  else if (!xorForceCms && !useXorCnf && pbBackendNeeded)
  {
    // Smart default: 'cvc' (Gauss-Jordan, ~27x faster) for LRA/theory inputs,
    // 'cms' (BV-PACT) for bit-vector or define-fun inputs where Gauss-Jordan
    // stalls or cannot run. An explicit --xor mode overrides this.
    useCvcXor = preferCvcXorDefault(filename);
  }
  // 'cms' / BV-PACT substitutes fresh 1-bit bit-vectors for the Boolean
  // projection variables so the XOR hashes reach the BV SAT backend.
  bool useBvPact = pbBackendNeeded && !useCvcXor && !useXorCnf
                   && !requestNativeXor;
  if (useCvcXor)
  {
    requestNativeXor = true;
    useBvPact = false;
    setenv("CVC5_XOR_GAUSS", "1", 1);
    // Enable Gauss-Jordan *propagation* (not just conflict detection) by
    // default: it prunes the search like CMS/ApproxMC and is the whole point
    // (~27x on case89). Set CVC5_XOR_GAUSS_NOPROP=1 to fall back to
    // conflict-only for A/B comparison.
    if (!getenv("CVC5_XOR_GAUSS_NOPROP"))
    {
      setenv("CVC5_XOR_GAUSS_PROP", "1", 1);
    }
  }
  if (useXorCnf)
  {
    // 'blast': TTC_XOR_CNF makes rebuildCountSolver force CaDiCaL with native
    // XOR disabled, so the parity is Tseitin-expanded to CNF and solved as
    // plain clauses -- no Gauss-Jordan anywhere.
    requestNativeXor = false;
    useBvPact = false;
    setenv("TTC_XOR_CNF", "1", 1);
  }
  bool useToCNF = vm.count("tocnf") > 0;
  bool useBvCnf = vm.count("bvcnf") > 0 || autoBvCnf;
  bool writeFB = vm.count("FB") > 0;
  bool useMono = vm.count("mono") > 0 || vm.count("mono-true") > 0;
  bool monoTrue = vm.count("mono-true") > 0;
#if defined(TTC_ENABLE_DDNNF)
  bool useResidualSimplifier = vm.count("no-res-simp") == 0;
  bool useD4 = vm.count("d4") > 0;
  bool fullDecompose = vm.count("full-decomp") > 0;
#else
  bool useResidualSimplifier = true;
  bool useD4 = false;
  bool fullDecompose = false;
#endif
  if (!ttc::ddnnfEnabled() && !useToCNF && !useBvCnf)
  {
    useApprox = true;
  }
  if (useProjectionBoost)
  {
    useApprox = true;
  }
  if (useProjectionEnumerate)
  {
    useApprox = false;
  }
  if (countUnsatAssignments && !useApprox)
  {
    std::cerr << "Error: --unsat is only supported with approximate counting (-a)"
              << std::endl;
    return 1;
  }
  std::string cnfFile;
  if (useToCNF) {
    if (!vm.count("output")) {
      std::cerr << "Error: --tocnf requires -o <file>" << std::endl;
      return 1;
    }
    cnfFile = vm["output"].as<std::string>();
  }
  std::string fbFile;
  if (writeFB) {
    fbFile = vm["FB"].as<std::string>();
  }
  bool doArjun = useArjun || useToCNF;
  bool useAssumptions = vm.count("assumptions") > 0;
  bool useBitset = vm.count("bitset-vars") > 0;
#if defined(TTC_ENABLE_DDNNF)
  bool useSmtCache = vm.count("smt-cache") > 0;
  bool noDecompose = vm.count("no-decompose") > 0;
  if (fullDecompose && noDecompose)
  {
    std::cerr << "Error: --full-decomp cannot be combined with --no-decompose"
              << std::endl;
    return 1;
  }
#else
  bool useSmtCache = false;
  bool noDecompose = false;
#endif
  int propAt = vm["propat"].as<int>();
  if (propAt < 1)
    propAt = 1;
  std::string cacheStr = vm["cache"].as<std::string>();
  CacheMode cacheMode = CacheMode::Syntactic;
  if (cacheStr == "hash") {
    cacheMode = CacheMode::Hash;
  } else if (cacheStr == "sem") {
    cacheMode = CacheMode::Semantic;
  } else if (cacheStr == "canon") {
    cacheMode = CacheMode::Canonical;
  } else if (cacheStr == "bool") {
    cacheMode = CacheMode::Bool;
  } else if (cacheStr != "synt") {
    std::cerr << "Error: Invalid cache mode '" << cacheStr << "'" << std::endl;
    return 1;
  }
  int verbosity = vm["verbose"].as<int>();
  if (verbosity < 0)
    verbosity = 0;
  if (verbosity > 10)
    verbosity = 10;
  Log.setVerbosity(verbosity);

  Trace.clear();
  if (vm.count("trace")) {
    auto traces = vm["trace"].as<std::vector<std::string>>();
    Trace.setEnabled(traces);
    if (!traces.empty()) {
      std::ostringstream oss;
      for (size_t i = 0; i < traces.size(); ++i) {
        if (i) oss << ", ";
        oss << traces[i];
      }
      Log(2) << "Enabled traces: " << oss.str() << std::endl;
    }
  }
  g_xorTraceEnabled = Trace.isEnabled("xor");

  double verbosityTime = vm["verbosity-time"].as<double>();
  if (verbosityTime <= 0)
    verbosityTime = 5.0;
  progress_interval = verbosityTime;

  Log(2) << "Opening input file: " << filename << std::endl;
  std::ifstream inFile(filename);
  if (!inFile) {
    std::cerr << "Error: Cannot open file " << filename << std::endl;
    return 1;
  }

  std::stringstream buffer;
  buffer << inFile.rdbuf();
  std::string smtFormula = buffer.str();
  if (useApprox)
  {
    smtFormula = ensureLogicIncludesBv(std::move(smtFormula));
  }
  if (useUnsatQuant)
  {
    smtFormula = overrideLogicForQuantifiers(std::move(smtFormula));
  }
  Log(2) << "Read input file" << std::endl;
  Log(3) << "File size: " << smtFormula.size() << " bytes" << std::endl;

  TTCParser parser;
  cvc5::Solver& mainSolver = parser.solver();
  if (std::getenv("TTC_PACT_STATS") != nullptr)
  {
    try { mainSolver.setOption("stats-internal", "true"); }
    catch (const cvc5::CVC5ApiException&) {}
  }
  try
  {
    mainSolver.setOption("bv-sat-solver", "cryptominisat");
  }
  catch (const cvc5::CVC5ApiException& ex)
  {
    Log(1) << "Warning: unable to select CryptoMiniSat backend: " << ex.getMessage()
           << std::endl;
  }
  if (useBvPact)
  {
    // Keep the fresh 1-bit projection bit-vectors as genuine bit-vectors so the
    // XOR hashes are solved by the BV SAT backend; otherwise cvc5 folds them
    // back into Booleans (bv-to-bool) and the core SAT solver handles them.
    try
    {
      mainSolver.setOption("bv-to-bool", "false");
    }
    catch (const cvc5::CVC5ApiException& ex)
    {
      Log(1) << "Warning: unable to disable bv-to-bool: " << ex.getMessage()
             << std::endl;
    }
  }


  g_hasCadicalXorSupport = detectCadicalXorSupport(mainSolver);
  if (requestNativeXor && !g_hasCadicalXorSupport)
  {
    std::cerr << "Warning: --xor requested but cvc5 lacks CaDiCaL native XOR support"
              << "; continuing without native XOR" << std::endl;
  }
  g_useNativeXor = requestNativeXor && g_hasCadicalXorSupport;
  if (g_hasCadicalXorSupport)
  {
    try
    {
      applyNativeXor(mainSolver, g_useNativeXor);
    }
    catch (const cvc5::CVC5ApiException& ex)
    {
      std::cerr << "Warning: unable to configure CaDiCaL native XOR support: "
                << ex.getMessage() << "; continuing without native XOR"
                << std::endl;
      g_useNativeXor = false;
      try
      {
        applyNativeXor(mainSolver, false);
      }
      catch (const cvc5::CVC5ApiException&)
      {
      }
    }
  }
  // mainSolver.setXorAssertionVerbose(true);
  if (useBvPact && g_hasCadicalXorSupport)
  {
    // BV-PACT routes the parity (XOR) hashes through the bit-vector SAT backend
    // (CryptoMiniSat).  Enabling native XOR makes cvc5's CNF stream hand each
    // bit-blasted XOR to the SAT solver as a single native XOR clause instead of
    // a Tseitin clause expansion.  This is independent of the (broken) ttc-level
    // hash_constraint native path, so we keep d_counter's useNativeXor disabled.
    try
    {
      applyNativeXor(mainSolver, true);
    }
    catch (const cvc5::CVC5ApiException& ex)
    {
      Log(1) << "Warning: unable to enable native XOR for BV-PACT: "
             << ex.getMessage() << std::endl;
    }
  }
  std::string solverInfo = "solver: cvc5";
  if (g_hasCadicalXorSupport)
  {
    solverInfo += " with cadical-xor";
  }
  set_banner_solver_info(std::move(solverInfo));
  Log(2) << "Parsing formula" << std::endl;
  double parseStart = Log.elapsed();
  parser.parseFormula(smtFormula);
  double parseEnd = Log.elapsed();
  Profile.addParse(parseEnd - parseStart);
  Log(3) << "Parsing complete" << std::endl;

  // When no projection variables were declared explicitly, default to projected
  // counting over the Boolean and bit-vector variables of the formula. This
  // mirrors cvc5's sampling-set behaviour for benchmarks that mix
  // Booleans/bit-vectors with other theories.
  if (!parser.hasExplicitProjectionVars() && parser.numProjVars() == 0)
  {
    parser.promoteBooleanAndBvToProjection();
  }

  // Resolve the deferred automatic engine choice now that the projection
  // variables (and their sorts) are known. --bvcnf (rule 1) was already
  // committed pre-parse.
  if (autoDispatch && autoDeferred)
  {
    const std::vector<cvc5::Term>& pv = parser.projectionVars();
    if (!pv.empty() && projVarsAreBvOrBool(pv))
    {
      // Rule 2: BV/Bool projection variables -> projection counting (the
      // pre-parse PB backend was already configured).
      useProjectionBoost = true;
      useApprox = true;
    }
    else
    {
      // Rule 3 candidate: pure LRA with no projection variables -> volume.
      // Mirror the isLraInput detection used below; anything else (e.g. a
      // bare QF_LIA, or non-BV/Bool projection variables) is unsupported.
      const bool lra =
          (!autoLogic.empty() && autoLogic.find("LRA") != std::string::npos)
          || (parser.numRealVars() > 0 && parser.numIntVars() == 0
              && parser.numBvVars() == 0);
      const bool volumeOk = lra && pv.empty() && !parser.formula().isNull();
      if (!volumeOk)
      {
        printUnsupportedEngine(autoLogic, parser.numProjVars());
        return 1;
      }
      // else: leave useProjectionBoost false; the isLraInput path runs volume.
    }
  }

  if (parser.hasWeights() && !useProjectionEnumerate && !useProjectionBoost)
  {
    std::cerr << "Error: declare-weight is currently supported only with --enum"
              << " or --pact" << std::endl;
    return 1;
  }
  if (parser.hasWeights() && countUnsatAssignments)
  {
    std::cerr << "Error: --unsat is not supported with declare-weight"
              << std::endl;
    return 1;
  }

  // Eager bit-blasting BV model counter: write a model-preserving CNF and run
  // ApproxMC over the bits of the bit-vector variables.
  if (useBvCnf)
  {
    const std::string bvCnfPath = vm.count("bvcnf")
                                      ? vm["bvcnf"].as<std::string>()
                                      : makeAutoCnfPath(filename);

    // Sampling set: every declared bit-vector variable, including those that
    // do not occur in any assertion (counted as free bits), so the projected
    // count matches the bit-vector model count of the formula.
    std::vector<cvc5::Term> bvVars;
    std::unordered_set<cvc5::Term> seenVars;
    for (const cvc5::Term& v : parser.declaredVariables())
    {
      if (v.getSort().isBitVector() && seenVars.insert(v).second)
      {
        bvVars.push_back(v);
      }
    }

    try
    {
      double bbStart = Log.elapsed();
      ttc::eager::BvCnfResult cnf = ttc::eager::writeBvCnf(
          parser.solver(), parser.assertions(), bvVars, bvCnfPath);
      double bbEnd = Log.elapsed();

      std::cout << "c [ttc->bvcnf] wrote model-preserving CNF to '" << cnf.path
                << "'" << std::endl;
      std::cout << "c [ttc->bvcnf] variables: " << cnf.numVars
                << " clauses: " << cnf.numClauses
                << " sampling: " << cnf.numSamplingVars << std::endl;
      std::cout << "c [ttc->bvcnf] bit-blasted in " << std::fixed
                << std::setprecision(2) << (bbEnd - bbStart) << " seconds"
                << std::endl;

      std::optional<std::string> count =
          ttc::eager::runApproxMc(bvCnfPath, approxSeed, verbosity);
      if (count.has_value())
      {
        std::cout << "s mc " << *count << std::endl;
        return 0;
      }
      std::cerr << "Error: ApproxMC did not produce a count" << std::endl;
      return 1;
    }
    catch (const std::exception& ex)
    {
      std::cerr << "Error: --bvcnf failed: " << ex.what() << std::endl;
      return 1;
    }
  }

  bool isLraInput = false;
  std::string logicName;
  try
  {
    logicName = mainSolver.getInfo("logic");
  }
  catch (const cvc5::CVC5ApiException&)
  {
    logicName.clear();
  }
  std::string logicUpper = logicName;
  std::transform(logicUpper.begin(), logicUpper.end(), logicUpper.begin(),
                 [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
  if (!logicUpper.empty() && logicUpper.find("LRA") != std::string::npos)
  {
    isLraInput = true;
  }
  else if (parser.numRealVars() > 0 && parser.numIntVars() == 0 &&
           parser.numBvVars() == 0)
  {
    isLraInput = true;
  }
  if (isLraInput && parser.formula().isNull())
  {
    isLraInput = false;
  }
  // Projection-based counting modes (--PB / --PE) count over the Boolean /
  // bit-vector support of the formula, not the LRA volume. Even for pure-LRA
  // logics they must run projected counting (pact) rather than fall into the
  // volume-computation path below.
  if (useProjectionBoost || useProjectionEnumerate)
  {
    isLraInput = false;
  }
  if (isLraInput)
  {
    try
    {
      print_banner();

      print_section("parsing input");
      std::cout << "c reading SMT2 file from '" << filename << "'" << std::endl;
      std::cout << "c parsed formula in " << std::fixed << std::setprecision(2)
                << (parseEnd - parseStart) << " seconds" << std::endl;
      std::cout << "c" << std::endl;

      double dnfStart = Log.elapsed();
      parser.computePolytopes();
      double dnfEnd = Log.elapsed();
      Profile.addDnfization(dnfEnd - dnfStart);

      const auto& polytopes = parser.polytopes();
      auto realVars = parser.realVariables();

      print_section("options");
      std::cout << "c counting: volume" << std::endl;
      std::cout << "c detected pure LRA: yes" << std::endl;
      std::cout << "c real variables: " << realVars.size() << std::endl;
      std::cout << "c polytopes: " << polytopes.size() << std::endl;
      std::cout << "c" << std::endl;

      print_section("volume computation");
      std::cout << "c polytope        volume   samples_deleted   total_samples"
                << std::endl;
      auto rowPrinter = [](const ttc::VolumeComputationRow& row) {
        std::cout << "c " << std::setw(8) << row.index << "  " << std::setw(14)
                  << std::fixed << std::setprecision(6) << row.volume << "  "
                  << std::setw(15) << row.samplesDeleted << "  "
                  << std::setw(14) << row.totalSamples << std::endl;
      };

      double volumeStart = Log.elapsed();
      auto volumeResult =
          ttc::computeLraVolume(polytopes, realVars, rowPrinter);
      double volumeEnd = Log.elapsed();
      Profile.addSearch(volumeEnd - volumeStart);
      Profile.addPolytopeVolume(volumeResult.volumeComputationTime);
      Profile.addSampling(volumeResult.samplingTime);
      Profile.setSampleStats(volumeResult.totalSamplesGenerated,
                             volumeResult.totalSamplesDeleted);
      std::cout.unsetf(std::ios::floatfield);
      std::cout.precision(6);
      std::cout << "c" << std::endl;

      print_section("result");
      std::cout << "s volume " << volumeResult.volumeEstimate << std::endl;
      std::cout << "c" << std::endl;

      double totalTime = Log.elapsed();
      Profile.print(totalTime);
      print_resources(totalTime);
      print_section("shutting down");
      std::cout << "c exit 0" << std::endl;
      return 0;
    }
    catch (const std::exception& ex)
    {
      std::cerr << "Error: volume computation failed: " << ex.what() << std::endl;
      return 1;
    }
  }

  auto stpRes = parser.checkSatWithSTP(smtFormula);
  if (stpRes.has_value()) {
    if (!stpRes.value()) {
      std::cout << "s unsat" << std::endl;
      return 0;
    } else {
      Log(2) << "STP satisfiable" << std::endl;
    }
  }
  if (doArjun) {
    double preStart = Log.elapsed();
    parser.minimizeProjectionSet();
    double preEnd = Log.elapsed();
    Profile.addPreprocess(preEnd - preStart);
  }

  if (writeFB) {
    ToCNF converter(parser.solver(), parser.assertions());
    auto cnf = converter.build();
    std::ofstream out(fbFile);
    out << "p cnf " << cnf.varCount << " " << cnf.clauses.size() << "\n";
    out << "c p show";
    for (const auto &t : parser.projectionVars()) {
      auto it = cnf.termToIdx.find(t);
      if (it != cnf.termToIdx.end()) {
        out << ' ' << it->second;
      }
    }
    out << " 0\n";
    for (const auto &cl : cnf.clauses) {
      for (int lit : cl)
        out << lit << ' ';
      out << "0\n";
    }
    out << "c mapping\n";
    for (int i = 1; i <= cnf.varCount; ++i) {
      out << "c " << i;
      if (i < static_cast<int>(cnf.idxToTerm.size()) &&
          !cnf.idxToTerm[i].isNull()) {
        out << " " << cnf.idxToTerm[i].toString();
      }
      out << "\n";
    }
    return 0;
  }

  if (useProjectionEnumerate)
  {
    ProjectedEnumerator enumerator(parser.solver(), parser.projectionVars());
    if (verbosity == 0)
    {
      if (parser.hasWeights())
      {
        long double res = enumerator.countWeighted(parser.literalWeights());
        std::cout << "s wmc ";
        printWeightedCount(res);
        std::cout << std::endl;
      }
      else
      {
        std::uint64_t res = enumerator.count();
        std::cout << "s mc " << res << std::endl;
      }
      return 0;
    }

    print_banner();

    print_section("parsing input");
    std::cout << "c reading SMT2 file from '" << filename << "'" << std::endl;
    std::cout << "c parsed formula in " << std::fixed << std::setprecision(2)
              << (parseEnd - parseStart) << " seconds" << std::endl;
    std::cout << "c" << std::endl;

    print_section("options");
    std::cout << "c counting: "
              << (parser.hasWeights() ? "weighted projected enumeration"
                                      : "projected enumeration")
              << std::endl;
    std::cout << "c projection boost: "
              << (useProjectionBoost ? "yes" : "no") << std::endl;
    std::cout << "c xor: " << (g_useNativeXor ? "yes" : "no") << std::endl;
    std::cout << "c" << std::endl;

    print_section("preprocessing");
    if (doArjun)
    {
      std::cout << "c projection vars before: " << parser.projVarsBeforeMin()
                << std::endl;
      std::cout << "c projection vars after : " << parser.projVarsAfterMin()
                << std::endl;
    }
    else
    {
      std::cout << "c projection vars: " << parser.numProjVars()
                << std::endl;
    }
    std::cout << "c" << std::endl;

    print_section("projected enumeration");
    double countStart = Log.elapsed();
    std::uint64_t res = 0;
    long double weightedRes = 0.0L;
    if (parser.hasWeights())
    {
      weightedRes = enumerator.countWeighted(parser.literalWeights());
    }
    else
    {
      res = enumerator.count();
    }
    double countEnd = Log.elapsed();
    Profile.addSearch(countEnd - countStart);
    std::uint64_t smtCalls = enumerator.getSmtCallCount();
    if (parser.hasWeights())
    {
      std::cout << "c weighted count   : ";
      printWeightedCount(weightedRes);
      std::cout << std::endl;
    }
    else
    {
      std::cout << "c enumerated models: " << res << std::endl;
    }
    std::cout << "c smt calls       : " << smtCalls << std::endl;
    std::cout << "c elapsed time    : " << std::fixed << std::setprecision(2)
              << (countEnd - countStart) << " seconds" << std::endl;
    std::cout.unsetf(std::ios::floatfield);
    std::cout.precision(6);
    std::cout << "c" << std::endl;

    print_section("result");
    if (parser.hasWeights())
    {
      std::cout << "s wmc ";
      printWeightedCount(weightedRes);
      std::cout << std::endl;
    }
    else
    {
      std::cout << "s mc " << res << std::endl;
    }
    std::cout << "c" << std::endl;

    double totalTime = Log.elapsed();
    print_resources(totalTime);
    print_section("shutting down");
    std::cout << "c exit 0" << std::endl;
    return 0;
  }

  if (useApprox) {
    if (verbosity == 0) {
      if (parser.hasWeights() && useProjectionBoost)
      {
        try
        {
          WeightedPbResult weighted =
              runWeightedProjectionBoost(parser,
                                         approxSeed,
                                         g_useNativeXor,
                                         useBvPact,
                                         approxEpsilon,
                                         approxDelta);
          std::cout << "s wmc ";
          printWeightedCount(weighted.weightedCount);
          std::cout << std::endl;
          return 0;
        }
        catch (const std::exception& ex)
        {
          std::cerr << "Error: weighted --pact failed: " << ex.what()
                    << std::endl;
          return 1;
        }
      }
      Pact counter(parser.solver(), parser.projectionVars(), approxSeed, g_useNativeXor, useBvPact, approxEpsilon, approxDelta);
      std::uint64_t res = counter.count();
      cpp_int unsatTotal;
      bool haveUnsat = false;
      if (countUnsatAssignments) {
        if (useUnsatQuant) {
          auto quantRes = countUnsatAssignmentsQuantified(parser, approxSeed);
          if (!quantRes.has_value()) {
            std::cerr << "Error: Unable to compute unsat count with quantifiers"
                      << std::endl;
            return 1;
          }
          unsatTotal = cpp_int(*quantRes);
          haveUnsat = true;
        }
        if (!haveUnsat) {
          auto totalAssignments = computeTotalAssignments(parser.projectionVars());
          if (!totalAssignments.has_value()) {
            std::cerr << "Error: --unsat currently supports only Boolean or"
                      << " bit-vector projection variables" << std::endl;
            return 1;
          }
          cpp_int satCount = cpp_int(res);
          unsatTotal = clampNonNegative(*totalAssignments - satCount);
          haveUnsat = true;
        }
      }
      std::cout << "s mc " << res << std::endl;
      if (countUnsatAssignments) {
        std::cout << "s unsat-mc " << unsatTotal << std::endl;
      }
      return 0;
    }
    markPactRunningInBanner();
    print_banner();

    print_section("parsing input");
    std::cout << "c reading SMT2 file from '" << filename << "'" << std::endl;
    std::cout << "c parsed formula in " << std::fixed << std::setprecision(2)
              << (parseEnd - parseStart) << " seconds" << std::endl;
    std::cout << "c" << std::endl;

    print_section("options");
    std::cout << "c counting: "
              << (parser.hasWeights() && useProjectionBoost
                      ? "weighted approximate"
                      : "approximate")
              << std::endl;
    std::cout << "c bv-pact: " << (useBvPact ? "yes (XOR hashes to BV SAT backend)"
                                             : "no (Boolean hashing backup)")
              << std::endl;
    std::cout << "c contract: " << (noContract ? "no" : "yes") << std::endl;
    std::cout << "c arjun: " << (doArjun ? "yes" : "no") << std::endl;
    std::cout << "c xor: " << (g_useNativeXor ? "yes" : "no") << std::endl;
    std::cout << "c unsat assignments: "
              << (countUnsatAssignments
                      ? (useUnsatQuant ? "yes (quantified)"
                                       : "yes (complement)")
                      : "no")
              << std::endl;
    std::cout << "c" << std::endl;

    print_section("preprocessing");
    if (doArjun) {
      std::cout << "c projection vars before: " << parser.projVarsBeforeMin()
                << std::endl;
      std::cout << "c projection vars after : " << parser.projVarsAfterMin()
                << std::endl;
    } else {
      std::cout << "c projection vars: " << parser.numProjVars() << std::endl;
    }
    std::cout << "c" << std::endl;

    print_section(parser.hasWeights() && useProjectionBoost
                      ? "weighted approximate counting"
                      : "approximate counting");
    // Column widths must match the data rows emitted by Pact::count so the
    // headings line up over their values.
    {
      std::ostringstream hdr;
      hdr << "c " << std::setw(7) << "sec";
      hdr << ' ' << std::setw(3) << "rnd";
      hdr << ' ' << std::setw(4) << "hash";
      hdr << ' ' << std::setw(8) << "sat";
      hdr << ' ' << std::setw(9) << "reuse";
      hdr << ' ' << std::setw(5) << "next";
      hdr << ' ' << std::setw(13) << "count";
      std::cout << hdr.str() << std::endl;
    }
    double countStart = Log.elapsed();
    std::uint64_t res = 0;
    WeightedPbResult weighted;
    try
    {
      if (parser.hasWeights() && useProjectionBoost)
      {
        weighted = runWeightedProjectionBoost(parser,
                                             approxSeed,
                                             g_useNativeXor,
                                             useBvPact,
                                             approxEpsilon,
                                             approxDelta);
      }
      else
      {
        Pact counter(parser.solver(), parser.projectionVars(), approxSeed, g_useNativeXor, useBvPact, approxEpsilon, approxDelta);
        res = counter.count();
        weighted.smtCalls = counter.getSmtCallCount();
      }
    }
    catch (const std::exception& ex)
    {
      std::cerr << "Error: approximate counting failed: " << ex.what()
                << std::endl;
      return 1;
    }
    double countEnd = Log.elapsed();
    Profile.addSearch(countEnd - countStart);
    std::uint64_t totalSmtCalls = weighted.smtCalls;
    cpp_int unsatTotal;
    bool haveUnsat = false;
    if (countUnsatAssignments) {
      if (useUnsatQuant) {
        print_section("unsat counting (quantified)");
        double unsatStart = Log.elapsed();
        std::uint64_t quantSmtCalls = 0;
        auto quantRes =
            countUnsatAssignmentsQuantified(parser, quantSmtCalls, approxSeed);
        double unsatEnd = Log.elapsed();
        if (!quantRes.has_value()) {
          std::cerr << "Error: Unable to compute unsat count with quantifiers"
                    << std::endl;
          return 1;
        }
        Profile.addSearch(unsatEnd - unsatStart);
        unsatTotal = cpp_int(*quantRes);
        haveUnsat = true;
        std::cout << "c" << std::endl;
        totalSmtCalls += quantSmtCalls;
      }
      if (!haveUnsat) {
        auto totalAssignments = computeTotalAssignments(parser.projectionVars());
        if (!totalAssignments.has_value()) {
          std::cerr << "Error: --unsat currently supports only Boolean or"
                    << " bit-vector projection variables" << std::endl;
          return 1;
        }
        cpp_int satCount = cpp_int(res);
        unsatTotal = clampNonNegative(*totalAssignments - satCount);
        haveUnsat = true;
      }
    }
    print_section("result");
    if (parser.hasWeights() && useProjectionBoost)
    {
      std::cout << "s wmc ";
      printWeightedCount(weighted.weightedCount);
      std::cout << std::endl;
    }
    else
    {
      std::cout << "s mc " << res << std::endl;
    }
    if (countUnsatAssignments) {
      std::cout << "s unsat-mc " << unsatTotal << std::endl;
    }
    std::cout << "c" << std::endl;

    print_section("statistics");
    if (parser.hasWeights() && useProjectionBoost)
    {
      std::cout << "c weighted-to-unweighted original vars: "
                << weighted.originalVars << std::endl;
      std::cout << "c weighted-to-unweighted original clauses: "
                << weighted.originalClauses << std::endl;
      std::cout << "c weighted-to-unweighted converted vars: "
                << weighted.convertedVars << std::endl;
      std::cout << "c weighted-to-unweighted converted clauses: "
                << weighted.convertedClauses << std::endl;
      std::cout << "c weighted-to-unweighted sampling vars: "
                << weighted.samplingVars << std::endl;
      std::cout << "c weighted-to-unweighted divisor: 2^"
                << weighted.divisorPower << std::endl;
      std::cout << "c weighted-to-unweighted multiplier: ";
      printWeightedCount(weighted.multiplier);
      std::cout << std::endl;
      std::cout << "c unweighted count: "
                << weighted.unweightedCount << std::endl;
    }
    std::cout << "c smt calls: " << totalSmtCalls << std::endl;
    std::cout << "c" << std::endl;

    double totalTime = Log.elapsed();
    print_section("profiling");
    std::cout << "c process time taken by individual solving procedures"
              << std::endl;
    std::cout << "c (percentage relative to process time for solving)"
              << std::endl;
    std::cout << "c" << std::endl;

    auto printApproxLine = [&](double t, double pct, const char* name) {
      std::cout << "c " << std::setw(10) << std::fixed << std::setprecision(2)
                << t;
      std::cout << "   " << std::setw(6) << std::fixed << std::setprecision(2)
                << pct << "% " << name << std::endl;
    };

    double parseTime = Profile.getParse();
    double preprocessTime = Profile.getPreprocess();
    double treeDecompTime = Profile.getTreeDecomp();
    double searchTime = Profile.getSearch();
    double rewriteTime = Profile.getRewrite();
    double componentTime = Profile.getComponent();
    double decomposeTime = Profile.getDecompose();
    double solveTime =
        parseTime + preprocessTime + treeDecompTime + searchTime;
    if (solveTime <= 0.0) {
      solveTime = 1e-9;
    }

    printApproxLine(parseTime, (parseTime / solveTime) * 100.0, "parse");
    if (preprocessTime > 0.0) {
      printApproxLine(preprocessTime,
                      (preprocessTime / solveTime) * 100.0,
                      "preprocessing");
    }
    if (treeDecompTime > 0.0) {
      printApproxLine(treeDecompTime,
                      (treeDecompTime / solveTime) * 100.0,
                      "tree-decomp");
    }
    double satCallTime = searchTime - rewriteTime - componentTime - decomposeTime;
    if (satCallTime < 0.0) {
      satCallTime = 0.0;
    }
    printApproxLine(satCallTime,
                    (satCallTime / solveTime) * 100.0,
                    "sat call");
    if (componentTime > 0.0) {
      printApproxLine(componentTime,
                      (componentTime / solveTime) * 100.0,
                      "components");
    }
    if (decomposeTime > 0.0) {
      printApproxLine(decomposeTime,
                      (decomposeTime / solveTime) * 100.0,
                      "decompose");
    }
    if (rewriteTime > 0.0) {
      printApproxLine(rewriteTime,
                      (rewriteTime / solveTime) * 100.0,
                      "rewrite");
    }
    std::cout << "c =================================" << std::endl;
    printApproxLine(solveTime, (solveTime / totalTime) * 100.0, "solve");
    std::cout << "c" << std::endl;
    std::cout << "c last line shows process time for solving" << std::endl;
    std::cout << "c (percentage relative to total process time)" << std::endl;
    std::cout << "c" << std::endl;

    print_resources(totalTime);
    print_section("shutting down");
    std::cout << "c exit 0" << std::endl;
    return 0;
  }

  auto runToCnf = [&]() -> int {
    if (verbosity == 0) {
      ToCNF converter(parser.solver(), parser.assertions());
      converter.convert(cnfFile);
      return 0;
    }

    print_banner();

    print_section("parsing input");
    std::cout << "c reading SMT2 file from '" << filename << "'" << std::endl;
    std::cout << "c parsed formula in " << std::fixed << std::setprecision(2)
              << (parseEnd - parseStart) << " seconds" << std::endl;
    std::cout << "c" << std::endl;

    print_section("options");
    std::cout << "c propagation: " << (simplifyProp ? "simplify" : "rewrite")
              << std::endl;
    std::cout << "c cache: "
              << (cacheMode == CacheMode::Hash
                      ? "hash"
                      : (cacheMode == CacheMode::Semantic
                             ? "semantic"
                             : (cacheMode == CacheMode::Canonical
                                    ? "canon"
                                    : (cacheMode == CacheMode::Bool
                                           ? "bool"
                                           : "syntactic"))))
              << std::endl;
    std::cout << "c smt-cache: " << (useSmtCache ? "yes" : "no") << std::endl;
    std::cout << "c no-decompose: " << (noDecompose ? "yes" : "no")
              << std::endl;
    std::cout << "c full-decomp: " << (fullDecompose ? "yes" : "no")
              << std::endl;
    std::cout << "c verbosity: " << verbosity << std::endl;
    std::cout << "c verbosity-time: " << verbosityTime << std::endl;
    std::cout << "c contract: " << (noContract ? "no" : "yes") << std::endl;
    std::cout << "c netrel: " << (useNetrel ? "yes" : "no") << std::endl;
    std::cout << "c arjun: " << (doArjun ? "yes" : "no") << std::endl;
    std::cout << "c assumptions: " << (useAssumptions ? "yes" : "no")
              << std::endl;
    std::cout << "c bitset-vars: " << (useBitset ? "yes" : "no") << std::endl;
    std::cout << "c propat: " << propAt << std::endl;
    std::cout << "c residual-simplifier: "
              << (useResidualSimplifier ? "yes" : "no") << std::endl;
    std::cout << "c mono: " << (useMono ? (monoTrue ? "true" : "false") : "off")
              << std::endl;
    std::cout << "c" << std::endl;

    print_section("preprocessing");
    if (doArjun) {
      std::cout << "c projection vars before: " << parser.projVarsBeforeMin()
                << std::endl;
      std::cout << "c projection vars after : " << parser.projVarsAfterMin()
                << std::endl;
    } else {
      std::cout << "c projection vars: " << parser.numProjVars() << std::endl;
    }
    std::cout << "c" << std::endl;

    print_section("statistics");
    std::cout << "c Boolean variables: " << parser.numBoolVars() << std::endl;
    std::cout << "c Bitvector variables: " << parser.numBvVars() << std::endl;
    std::cout << "c Integer variables: " << parser.numIntVars() << std::endl;
    std::cout << "c Real variables: " << parser.numRealVars() << std::endl;
    std::cout << "c constraints: " << parser.numConstraints() << std::endl;
    std::cout << "c" << std::endl;

    parser.printTreeDecomposition(!noContract, useNetrel);

    print_section("converting", false);
    std::cout << "c #assert    sec    #smt_op    #bool_var    #cnf_clauses"
              << std::endl;
    ToCNF converter(parser.solver(), parser.assertions());
    auto stats = converter.convert(cnfFile);
    for (std::size_t i = 0; i < stats.size(); ++i) {
      const auto &st = stats[i];
      std::cout << "c " << std::setw(7) << (i + 1) << " " << std::setw(7)
                << std::fixed << std::setprecision(2) << st.timeSec << " "
                << std::setw(10) << st.smtOps << " " << std::setw(10)
                << st.boolVars << " " << std::setw(12) << st.clauses
                << std::endl;
    }
    std::cout << "c" << std::endl;

    print_section("result");
    std::cout << "c written CNF to " << cnfFile << std::endl;
    std::cout << "c" << std::endl;

    double totalTime = Log.elapsed();
    Log(3) << "Total time: " << totalTime << " s" << std::endl;
    Profile.print(totalTime);
    print_resources(totalTime);
    print_section("shutting down");
    std::cout << "c exit 0" << std::endl;
    return 0;
  };

#if defined(TTC_ENABLE_DDNNF)
  if (useToCNF)
  {
    return runToCnf();
  }
  else if (useD4)
  {
    ToCNF converter(parser.solver(), parser.assertions());
    auto cnf = converter.build();
    std::vector<int> projIds;
    for (const auto& t : parser.projectionVars())
    {
      auto it = cnf.termToIdx.find(t);
      if (it != cnf.termToIdx.end())
      {
        projIds.push_back(it->second);
      }
      else
      {
        auto strIt = cnf.varToIdx.find(t.toString());
        if (strIt != cnf.varToIdx.end())
        {
          projIds.push_back(strIt->second);
        }
      }
    }

    auto& parserTm = ttc::getTermBuilder(parser.solver());
    cvc5::Term formula =
        parser.assertions().size() == 1
            ? parser.assertions()[0]
            : parserTm.mkTerm(cvc5::Kind::AND, parser.assertions());
    auto orderedProj =
        computeProjVarOrder(formula, parser.projectionVars(), parser.solver(),
                            false, !noContract, useNetrel);
    std::vector<double> varOrder(cnf.varCount + 1,
                                 static_cast<double>(cnf.varCount + 1));
    for (std::size_t i = 0; i < orderedProj.size(); ++i)
    {
      auto it = cnf.termToIdx.find(orderedProj[i]);
      if (it != cnf.termToIdx.end())
      {
        varOrder[it->second] = static_cast<double>(i);
      }
      else
      {
        auto strIt = cnf.varToIdx.find(orderedProj[i].toString());
        if (strIt != cnf.varToIdx.end())
        {
          varOrder[strIt->second] = static_cast<double>(i);
        }
      }
    }

    if (verbosity == 0)
    {
      std::uint64_t count = d4Count(
          cnf.varCount, cnf.clauses, projIds, varOrder, parser.solver(),
          cnf.idxToTerm, useSmtCache, noDecompose, useResidualSimplifier,
          fullDecompose, useMono, monoTrue, cacheMode == CacheMode::Bool,
          nullptr);
      std::cout << "s mc " << count << std::endl;
      return 0;
    }

    print_banner();

    print_section("parsing input");
    std::cout << "c reading SMT2 file from '" << filename << "'" << std::endl;
    std::cout << "c parsed formula in " << std::fixed << std::setprecision(2)
              << (parseEnd - parseStart) << " seconds" << std::endl;
    std::cout << "c" << std::endl;

    print_section("options");
    std::cout << "c propagation: " << (simplifyProp ? "simplify" : "rewrite")
              << std::endl;
    std::cout << "c cache: "
              << (cacheMode == CacheMode::Hash
                      ? "hash"
                      : (cacheMode == CacheMode::Semantic
                             ? "semantic"
                             : (cacheMode == CacheMode::Canonical
                                    ? "canon"
                                    : (cacheMode == CacheMode::Bool
                                           ? "bool"
                                           : "syntactic"))))
              << std::endl;
    std::cout << "c smt-cache: " << (useSmtCache ? "yes" : "no") << std::endl;
    std::cout << "c no-decompose: " << (noDecompose ? "yes" : "no")
              << std::endl;
    std::cout << "c full-decomp: " << (fullDecompose ? "yes" : "no")
              << std::endl;
    std::cout << "c verbosity: " << verbosity << std::endl;
    std::cout << "c verbosity-time: " << verbosityTime << std::endl;
    std::cout << "c contract: " << (noContract ? "no" : "yes") << std::endl;
    std::cout << "c netrel: " << (useNetrel ? "yes" : "no") << std::endl;
    std::cout << "c arjun: " << (doArjun ? "yes" : "no") << std::endl;
    std::cout << "c assumptions: " << (useAssumptions ? "yes" : "no")
              << std::endl;
    std::cout << "c bitset-vars: " << (useBitset ? "yes" : "no") << std::endl;
    std::cout << "c propat: " << propAt << std::endl;
    std::cout << "c mono: " << (useMono ? (monoTrue ? "true" : "false") : "off")
              << std::endl;
    std::cout << "c" << std::endl;

    print_section("preprocessing");
    if (doArjun)
    {
      std::cout << "c projection vars before: " << parser.projVarsBeforeMin()
                << std::endl;
      std::cout << "c projection vars after : " << parser.projVarsAfterMin()
                << std::endl;
    }
    else
    {
      std::cout << "c projection vars: " << parser.numProjVars() << std::endl;
    }
    std::cout << "c" << std::endl;

    print_section("statistics");
    std::cout << "c Boolean variables: " << parser.numBoolVars() << std::endl;
    std::cout << "c Bitvector variables: " << parser.numBvVars() << std::endl;
    std::cout << "c Integer variables: " << parser.numIntVars() << std::endl;
    std::cout << "c Real variables: " << parser.numRealVars() << std::endl;
    std::cout << "c constraints: " << parser.numConstraints() << std::endl;
    std::cout << "c" << std::endl;

    parser.printTreeDecomposition(!noContract, useNetrel);

    print_section("d4 counting", false);
    std::cout << "c cnf variables: " << cnf.varCount
              << ", clauses: " << cnf.clauses.size()
              << ", projection vars: " << projIds.size() << std::endl;

    std::cout << "c" << std::endl;
    std::cout << "c   sec   cache_hit/miss  smtcall(unsat)  decisions   decomp(stopped)"
              << std::endl;

    Profile.resetSearchStats();
    if (cacheMode == CacheMode::Bool)
    {
      Profile.setCacheStats(0, 0);
    }
    double countStart = Log.elapsed();
    start_progress(countStart);
    D4Statistics d4Stats;
    std::uint64_t count = d4Count(
        cnf.varCount, cnf.clauses, projIds, varOrder, parser.solver(),
        cnf.idxToTerm, useSmtCache, noDecompose, useResidualSimplifier,
        fullDecompose, useMono, monoTrue, cacheMode == CacheMode::Bool,
        &d4Stats);
    double countEnd = Log.elapsed();
    Profile.addSearch(countEnd - countStart);
    if (cacheMode == CacheMode::Bool)
    {
      Profile.setCacheStats(d4Stats.cacheHits + d4Stats.cacheMisses,
                            d4Stats.cacheHits);
    }
    if (cacheMode == CacheMode::Bool)
    {
      double elapsed = countEnd - countStart;
      std::uint64_t smtCalls = Profile.getCheckSatCalls();
      std::uint64_t smtUnsat = Profile.getUnsatCheckSatCalls();
      std::uint64_t decisions = Profile.getDecisions();
      std::uint64_t decomps = Profile.getDecompositions();
      std::uint64_t stopped = Profile.getStoppedDecompositions();
      std::cout << "c " << std::fixed << std::setprecision(3) << elapsed
                << "  " << d4Stats.cacheHits << "/" << d4Stats.cacheMisses
                << "  " << smtCalls << "(" << smtUnsat << ")"
                << "  " << decisions
                << "  " << decomps << "(" << stopped << ")" << std::endl;
      std::cout.unsetf(std::ios::floatfield);
      std::cout.precision(6);
    }
    std::cout << "c" << std::endl;

    print_section("result");
    std::cout << "s mc " << count << std::endl;
    std::cout << "c" << std::endl;

    double totalTime = Log.elapsed();
    Log(3) << "Total time: " << totalTime << " s" << std::endl;
    Profile.print(totalTime);
    print_resources(totalTime);
    print_section("shutting down");
    std::cout << "c exit 0" << std::endl;
    return 0;
  }
  else if (verbosity == 0)
  {
    std::uint64_t count = parser.projectedModelCount(
        cacheMode, simplifyProp, !noContract, useNetrel, useAssumptions,
        useBitset, propAt, useMono, monoTrue);
    std::cout << "s mc " << count << std::endl;
    return 0;
  }

  print_banner();

  print_section("parsing input");
  std::cout << "c reading SMT2 file from '" << filename << "'" << std::endl;
  std::cout << "c parsed formula in " << std::fixed << std::setprecision(2)
            << (parseEnd - parseStart) << " seconds" << std::endl;
  std::cout << "c" << std::endl;

  print_section("options");
  std::cout << "c propagation: " << (simplifyProp ? "simplify" : "rewrite")
            << std::endl;
  std::cout << "c cache: "
            << (cacheMode == CacheMode::Hash
                    ? "hash"
                    : (cacheMode == CacheMode::Semantic
                           ? "semantic"
                           : (cacheMode == CacheMode::Canonical
                                  ? "canon"
                                  : (cacheMode == CacheMode::Bool
                                         ? "bool"
                                         : "syntactic"))))
            << std::endl;
  std::cout << "c smt-cache: " << (useSmtCache ? "yes" : "no") << std::endl;
  std::cout << "c no-decompose: " << (noDecompose ? "yes" : "no") << std::endl;
  std::cout << "c full-decomp: " << (fullDecompose ? "yes" : "no")
            << std::endl;
  std::cout << "c verbosity: " << verbosity << std::endl;
  std::cout << "c verbosity-time: " << verbosityTime << std::endl;
  std::cout << "c contract: " << (noContract ? "no" : "yes") << std::endl;
  std::cout << "c netrel: " << (useNetrel ? "yes" : "no") << std::endl;
  std::cout << "c arjun: " << (doArjun ? "yes" : "no") << std::endl;
  std::cout << "c assumptions: " << (useAssumptions ? "yes" : "no")
            << std::endl;
  std::cout << "c bitset-vars: " << (useBitset ? "yes" : "no") << std::endl;
  std::cout << "c propat: " << propAt << std::endl;
  std::cout << "c residual-simplifier: "
            << (useResidualSimplifier ? "yes" : "no") << std::endl;
  std::cout << "c mono: " << (useMono ? (monoTrue ? "true" : "false") : "off")
            << std::endl;
  std::cout << "c" << std::endl;

  print_section("preprocessing");
  if (doArjun)
  {
    std::cout << "c projection vars before: " << parser.projVarsBeforeMin()
              << std::endl;
    std::cout << "c projection vars after : " << parser.projVarsAfterMin()
              << std::endl;
  }
  else
  {
    std::cout << "c projection vars: " << parser.numProjVars() << std::endl;
  }
  std::cout << "c" << std::endl;

  print_section("statistics");
  std::cout << "c Boolean variables: " << parser.numBoolVars() << std::endl;
  std::cout << "c Bitvector variables: " << parser.numBvVars() << std::endl;
  std::cout << "c Integer variables: " << parser.numIntVars() << std::endl;
  std::cout << "c Real variables: " << parser.numRealVars() << std::endl;
  std::cout << "c constraints: " << parser.numConstraints() << std::endl;
  std::cout << "c" << std::endl;

  parser.printTreeDecomposition(!noContract, useNetrel);

  print_section("counting", false);
  std::cout
      << "c     sec      percentage     cache_size     lookups     cache_hit   "
         "decompositions   mono_clauses   check_sat   decisions"
      << std::endl;
  double countStart = Log.elapsed();
  start_progress(countStart);
  std::uint64_t count = parser.projectedModelCount(
      cacheMode, simplifyProp, !noContract, useNetrel, useAssumptions,
      useBitset, propAt, useMono, monoTrue);
  double countEnd = Log.elapsed();
  Profile.addSearch(countEnd - countStart);
  Log(3) << "Counting finished" << std::endl;
  std::cout << "c" << std::endl;

  print_section("result");
  std::cout << "s mc " << count << std::endl;
  std::cout << "c" << std::endl;

  double totalTime = Log.elapsed();
  Log(3) << "Total time: " << totalTime << " s" << std::endl;
  Profile.print(totalTime);

  print_resources(totalTime);

  print_section("shutting down");
  std::cout << "c exit 0" << std::endl;

  return 0;
#else
  if (useToCNF)
  {
    return runToCnf();
  }

  std::cerr << "Error: DDNNF/D4 counting is disabled in this build."
            << " Reconfigure with --ddnnf or use approximate counting." << std::endl;
  return 1;
#endif
}
