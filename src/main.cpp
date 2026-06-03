#include <algorithm>
#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <optional>
#include <utility>
#include <cctype>

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

}  // namespace

int main(int argc, char *argv[]) {
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "Display help information")
      ("rewrite-prop", "Use rewrite-based propagation instead of cvc5 simplify")
      ("no-contract",
       "Do not contract to projection variables for tree decomposition")
      ("netrel", "Apply netlist relevance reduction before contraction")
      ("arjun", "Minimize projection set before counting")
      ("approx,a",
       "Use approximate model counting")
      ("seed,s",
       po::value<std::uint64_t>()->default_value(42),
       "Seed for randomized approximate counting")
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
      ("unsat",
       "Count projection assignments that make the constraints unsatisfiable")
      ("unsat-q",
       "Use quantified reasoning when counting unsatisfied projection assignments")
      ("PB", "Enable projection-based approximate counting mode")
      ("no-bv-pact",
       "In --PB mode, keep Boolean projection variables (and route XOR hashes "
       "through the core SAT solver) instead of substituting fresh 1-bit "
       "bit-vector variables that send the XOR hashes to the BV SAT backend")
      ("PE", "Enumerate projected assignments exactly")
      ("xor",
       po::bool_switch()->default_value(false),
       "Enable CaDiCaL native XOR support when available")
      ("verbose,v", po::value<int>()->default_value(1),
       "Verbosity level (0-10)")
      ("verbosity-time,vt",
       po::value<double>()->default_value(5.0),
       "Print progress every K seconds")
      ("trace,t",
       po::value<std::vector<std::string>>()->composing()->multitoken(),
       "Enable tracing for selected components")
      ("input-file", po::value<std::string>(), "SMT2 file");

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
      else if (a == "-PB")
        a = "--PB";
      else if (a == "-PE")
        a = "--PE";
      args.push_back(a);
    }
    po::store(po::command_line_parser(args).options(desc).positional(p).run(),
              vm);
    po::notify(vm);
  } catch (const std::exception &ex) {
    std::cerr << "Error parsing command line: " << ex.what() << std::endl;
    return 1;
  }

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << "Usage: " << argv[0] << " [options] filename.smt2"
              << std::endl;
    std::cout << desc << std::endl;
    return 1;
  }

  std::string filename = vm["input-file"].as<std::string>();
  bool simplifyProp = vm.count("rewrite-prop") == 0;
  bool noContract = vm.count("no-contract") > 0;
  bool useNetrel = vm.count("netrel") > 0;
  bool useArjun = vm.count("arjun") > 0;
  bool useApprox = vm.count("approx") > 0;
  bool useProjectionBoost = vm.count("PB") > 0;
  bool useProjectionEnumerate = vm.count("PE") > 0;
  // BV-PACT is the default for --PB: substitute fresh 1-bit bit-vectors for the
  // Boolean projection variables so the XOR hashes reach the BV SAT backend.
  // --no-bv-pact restores the Boolean-hashing backup for cross-checking.
  bool useBvPact = useProjectionBoost && vm.count("no-bv-pact") == 0;
  std::uint64_t approxSeed = vm["seed"].as<std::uint64_t>();
  bool useUnsatQuant = vm.count("unsat-q") > 0;
  bool countUnsatAssignments = vm.count("unsat") > 0 || useUnsatQuant;
  bool requestNativeXor = vm["xor"].as<bool>();
  bool useToCNF = vm.count("tocnf") > 0;
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
  if (!ttc::ddnnfEnabled() && !useToCNF)
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

  // When no explicit projection variables (prefixed "proj_") were declared,
  // default to projected counting over the Boolean and bit-vector variables of
  // the formula. This mirrors cvc5's sampling-set behaviour for benchmarks that
  // mix Booleans/bit-vectors with other theories.
  if (parser.numProjVars() == 0)
  {
    parser.promoteBooleanAndBvToProjection();
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
      std::uint64_t res = enumerator.count();
      std::cout << "s mc " << res << std::endl;
      return 0;
    }

    print_banner();

    print_section("parsing input");
    std::cout << "c reading SMT2 file from '" << filename << "'" << std::endl;
    std::cout << "c parsed formula in " << std::fixed << std::setprecision(2)
              << (parseEnd - parseStart) << " seconds" << std::endl;
    std::cout << "c" << std::endl;

    print_section("options");
    std::cout << "c counting: projected enumeration" << std::endl;
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
    std::uint64_t res = enumerator.count();
    double countEnd = Log.elapsed();
    Profile.addSearch(countEnd - countStart);
    std::uint64_t smtCalls = enumerator.getSmtCallCount();
    std::cout << "c enumerated models: " << res << std::endl;
    std::cout << "c smt calls       : " << smtCalls << std::endl;
    std::cout << "c elapsed time    : " << std::fixed << std::setprecision(2)
              << (countEnd - countStart) << " seconds" << std::endl;
    std::cout.unsetf(std::ios::floatfield);
    std::cout.precision(6);
    std::cout << "c" << std::endl;

    print_section("result");
    std::cout << "s mc " << res << std::endl;
    std::cout << "c" << std::endl;

    double totalTime = Log.elapsed();
    print_resources(totalTime);
    print_section("shutting down");
    std::cout << "c exit 0" << std::endl;
    return 0;
  }

  if (useApprox) {
    if (verbosity == 0) {
      Pact counter(parser.solver(), parser.projectionVars(), approxSeed, g_useNativeXor, useBvPact);
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
    std::cout << "c counting: approximate" << std::endl;
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

    print_section("approximate counting");
    // Column widths must match the data rows emitted by Pact::count so the
    // headings line up over their values.
    {
      std::ostringstream hdr;
      hdr << "c " << std::setw(8) << "sec";
      hdr << "  " << std::setw(6) << "round";
      hdr << "  " << std::setw(6) << "hashes";
      hdr << "  " << std::setw(10) << "saturating";
      hdr << "  " << std::setw(9) << "nexthash";
      hdr << "  " << std::setw(16) << "count";
      std::cout << hdr.str() << std::endl;
    }
    Pact counter(parser.solver(), parser.projectionVars(), approxSeed, g_useNativeXor, useBvPact);
    double countStart = Log.elapsed();
    std::uint64_t res = counter.count();
    double countEnd = Log.elapsed();
    Profile.addSearch(countEnd - countStart);
    std::uint64_t totalSmtCalls = counter.getSmtCallCount();
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
    std::cout << "s mc " << res << std::endl;
    if (countUnsatAssignments) {
      std::cout << "s unsat-mc " << unsatTotal << std::endl;
    }
    std::cout << "c" << std::endl;

    print_section("statistics");
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
