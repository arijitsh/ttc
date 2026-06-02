#include "d4_interface.hpp"

#include <boost/multiprecision/gmp.hpp>
#include <iostream>
#include <sstream>

#include "d4/src/configurations/ConfigurationProjMcMethod.hpp"
#include "d4/src/methods/ProjMCMethod.hpp"
#include "d4/src/options/methods/OptionProjMcMethod.hpp"
#include "d4/src/options/branchingHeuristic/OptionBranchingHeuristic.hpp"
#include "d4/src/problem/cnf/ProblemManagerCnf.hpp"
#include "d4_smt.hpp"
#include "ddnnf/logger.hpp"

using D4Cache = ComponentCache<cvc5::Term, std::hash<cvc5::Term>>;

cvc5::Solver *d4_smt_solver = nullptr;
std::vector<cvc5::Term> d4_smt_idxToTerm;
std::unordered_set<int> d4_smt_projVars;
bool d4_use_mono = false;
bool d4_no_decompose = false;
bool d4_full_decompose = false;
D4Cache d4_smt_component_cache(true, [](const cvc5::Term &, const cvc5::Term &){ return false; });
D4Cache *d4_component_cache = nullptr;
bool d4_use_component_cache = false;
bool d4_use_residual_simplifier = true;
std::unordered_map<const void*,
                   std::unordered_map<std::size_t, boost::multiprecision::mpz_int>>
    d4_component_smt_values;

std::uint64_t d4Count(int numVars, const std::vector<std::vector<int>> &cnf,
                      const std::vector<int> &projVars,
                      const std::vector<double> &varOrder, cvc5::Solver &solver,
                      const std::vector<cvc5::Term> &idxToTerm, bool useSmtCache,
                      bool noDecompose, bool useResidualSimplifier,
                      bool fullDecompose, bool useMono, bool monoTrue,
                      bool useComponentCache, D4Statistics *stats) {
  using namespace d4;
  using boost::multiprecision::mpf_float;
  using boost::multiprecision::mpz_int;

  if (stats) {
    stats->cacheHits = 0;
    stats->cacheMisses = 0;
  }

  d4_smt_solver = &solver;
  d4_no_decompose = noDecompose;
  d4_use_residual_simplifier = useResidualSimplifier;
  d4_full_decompose = fullDecompose;
  d4_use_component_cache = useComponentCache;
  // Copy the index-to-term mapping directly so that MiniSAT variable numbers
  // align with their corresponding SMT terms.  Both MiniSAT and the CNF
  // abstraction use 1-based variable identifiers, leaving index 0 unused.
  // Hence we keep the dummy 0th element rather than shifting the vector, which
  // previously caused off-by-one errors when synchronising assignments.
  d4_smt_idxToTerm = idxToTerm;
  if (Trace.isEnabled("mapping")) {
    Trace("mapping") << "c Boolean abstraction mapping (var -> term)" << std::endl;
    for (std::size_t idx = 1; idx < d4_smt_idxToTerm.size(); ++idx) {
      const cvc5::Term &term = d4_smt_idxToTerm[idx];
      if (term.isNull()) continue;
      Trace("mapping") << "c   " << idx << " -> " << term << std::endl;
    }
    std::vector<cvc5::Term> assertions = solver.getAssertions();
    if (!assertions.empty()) {
      Trace("mapping") << "c SMT assertions:" << std::endl;
      for (std::size_t i = 0; i < assertions.size(); ++i) {
        Trace("mapping") << "c   [" << i << "] " << assertions[i] << std::endl;
      }
    }
  }
  d4_smt_projVars.clear();
  // Store projection variables using the same indexing as the CNF/solver.
  // The earlier implementation subtracted 1 from each identifier, resulting in
  // checks for a non-existent variable 0 and preventing SMT consistency checks
  // from ever being triggered.
  for (int v : projVars)
    d4_smt_projVars.insert(v);

  std::vector<mpf_float> weightLit((numVars + 1) * 2, mpf_float(1));
  std::vector<mpf_float> weightVar(numVars + 1, mpf_float(1));
  std::vector<Var> selected;
  selected.reserve(projVars.size());
  for (int v : projVars) {
    selected.push_back(v);
    weightVar[v] = mpf_float(2);
  }

  auto *pm = new ProblemManagerCnf(numVars, weightLit, weightVar, selected);

  std::vector<std::vector<Lit>> d4Clauses;
  d4Clauses.reserve(cnf.size());
  for (const auto &cl : cnf) {
    std::vector<Lit> dcl;
    dcl.reserve(cl.size());
    for (int lit : cl) {
      bool sign = lit < 0;
      int v = sign ? -lit : lit;
      dcl.push_back(Lit::makeLit(v, sign));
    }
    d4Clauses.push_back(std::move(dcl));
  }
  pm->setClauses(d4Clauses);

  ConfigurationProjMcMethod config;
  config.counter.branchingHeuristic.configurationPartialOrderHeuristic
      .partialOrderMethod = d4::PARTIAL_ORDER_GIVEN;
  config.counter.branchingHeuristic.configurationPartialOrderHeuristic
      .givenOrder = varOrder;
  config.counter.branchingHeuristic.configurationPartialOrderHeuristic
      .scaleFactor = static_cast<double>(numVars + 1);
  config.counter.verbosity = false;
  if (useMono) {
    config.counter.branchingHeuristic.phaseHeuristicType =
        monoTrue ? d4::PHASE_TRUE : d4::PHASE_FALSE;
  }
  // Enable the native component cache within the DPLL-style counter when
  // requested.  The outer ProjMC cache remains disabled because the SMT-aware
  // integration only supports caching at the counter level.
  config.cache.isActivated = false;
  config.counter.cache.isActivated = useComponentCache;
  if (noDecompose) {
    config.counter.decompose = false;
  }
  OptionProjMcMethod opts(config);
  std::ostringstream d4Log;
  d4_use_mono = useMono;
  if (useSmtCache)
    d4_component_cache = &d4_smt_component_cache;
  else
    d4_component_cache = nullptr;
  ProjMCMethod<mpz_int> method(opts, pm, d4Log);
  mpz_int res = method.run();

  if (stats) {
    stats->cacheHits = method.cachePositiveHits();
    stats->cacheMisses = method.cacheNegativeHits();
  }
  // The SMT terms stored in the component cache reference the solver's term
  // manager.  Clear the cache before the solver is destroyed to avoid dangling
  // references during static destruction at program shutdown.
  if (useSmtCache) {
    d4_smt_component_cache.clear();
  }
  d4_component_smt_values.clear();

  d4_smt_solver = nullptr;
  d4_smt_idxToTerm.clear();
  d4_smt_projVars.clear();
  d4_use_mono = false;
  d4_no_decompose = false;
  d4_use_residual_simplifier = true;
  d4_full_decompose = false;
  d4_use_component_cache = false;
  d4_component_cache = nullptr;
  return res.convert_to<std::uint64_t>();
}
