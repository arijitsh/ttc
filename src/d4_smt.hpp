#pragma once
#include <cvc5/cvc5.h>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <cstddef>
#include <boost/multiprecision/gmp.hpp>
#include "ddnnf/cache.hpp"

// Global hooks for synchronising D4's Boolean search with the original SMT formula.
extern cvc5::Solver* d4_smt_solver;
extern std::vector<cvc5::Term> d4_smt_idxToTerm;
extern std::unordered_set<int> d4_smt_projVars;
// Enable monotone clause learning inside the D4 backend when set to true.
extern bool d4_use_mono;
// When true, the D4 backend skips component-aware bookkeeping to mimic
// the behaviour of the legacy implementation without decomposition.
extern bool d4_no_decompose;
// When true, decomposition ignores SMT-level dependencies.
extern bool d4_full_decompose;
// Disable the residual SMT simplifier when false.
extern bool d4_use_residual_simplifier;
// SMT-aware component cache used by the D4 backend when enabled.
extern ComponentCache<cvc5::Term, std::hash<cvc5::Term>>* d4_component_cache;
// Flag indicating whether D4's native component cache is active for the
// current run.
extern bool d4_use_component_cache;
// Mapping from cached component representations to SMT-aware counts.  When the
// SMT context differs, components reuse the cached structure but store
// disambiguated counts keyed by a hash of the SMT assignments.
extern std::unordered_map<const void*,
                          std::unordered_map<std::size_t, boost::multiprecision::mpz_int>>
    d4_component_smt_values;
