#ifndef PROJ_DDNNF_HPP
#define PROJ_DDNNF_HPP

#include <cvc5/cvc5.h>
#include <cstdint>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <string>
#include "bitset.hpp"
#include "cache.hpp"
#include "cache_mode.hpp"

// A simple DDNNF-like model counter that decomposes the formula
// when possible and branches on projection variables when decomposition is not possible.
class ProjDDNNF {
public:
    ProjDDNNF(cvc5::Solver& solver,
                const cvc5::Term& formula,
                const std::vector<cvc5::Term>& projVars,
                CacheMode cacheMode,
                bool simplifyProp,
                bool contract,
                bool netrel,
                bool useAssumptions,
                bool useBitset,
                int propAt,
                bool useMono,
                bool monoTrue);

    std::uint64_t count();

private:
    std::uint64_t countRec(const cvc5::Term& formula,
                           const std::vector<cvc5::Term>& projVars,
                           std::size_t depth,
                           std::vector<std::pair<std::size_t, bool>>& assigns);
    cvc5::Term assignAndRewrite(const cvc5::Term& formula,
                                const cvc5::Term& var,
                                bool value,
                                std::size_t depth);
    void collectVars(const cvc5::Term& term,
                     std::unordered_set<std::string>& vars);
    void collectVars(const cvc5::Term& term, sspp::Bitset& vars);
    void addMonotoneClause();

    cvc5::Solver& d_solver;
    cvc5::Term d_formula;
    std::vector<cvc5::Term> d_projVars;
    double d_totalNodes = 0;
    double d_nodesVisited = 0;
    bool d_simplifyProp;
    bool d_useAssumptions;
    bool d_useBitset;
    CacheMode d_cacheMode;
    std::vector<cvc5::Term> d_assumptions;
    int d_propatK;
    bool d_useMono;
    bool d_monoTrue;
    std::vector<cvc5::Term> d_falseStack;
    std::unordered_map<cvc5::Term, std::size_t, std::hash<cvc5::Term>> d_varToId;
    std::size_t d_varCount = 0;
    ComponentCache<DDNNFCacheKey, DDNNFCacheKeyHash> d_cache;
    std::vector<std::size_t> extractConstraintIndices(const cvc5::Term& formula);
};

#endif // PROJ_DDNNF_HPP
