#include "proj_ddnnf.hpp"

#include <functional>
#include <numeric>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <array>

#include "features.hpp"
#include "logger.hpp"
#include "profiler.hpp"
#include "rewrite.hpp"
#include "var_order.hpp"

ProjDDNNF::ProjDDNNF(cvc5::Solver& solver,
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
                       bool monoTrue)
    : d_solver(solver),
      d_formula(formula),
      d_projVars(computeProjVarOrder(formula, projVars, solver, false, contract, netrel)),
      d_simplifyProp(simplifyProp),
      d_useAssumptions(useAssumptions),
      d_useBitset(useBitset),
      d_cacheMode(cacheMode),
      d_propatK(propAt),
      d_useMono(useMono),
      d_monoTrue(monoTrue),
      d_cache(cacheMode == CacheMode::Hash || cacheMode == CacheMode::Bool,
              [&solver, cacheMode](const DDNNFCacheKey& a, const DDNNFCacheKey& b) {
                  if (a.vars != b.vars)
                  {
                      return false;
                  }
                  if (a.formula == b.formula)
                  {
                      return true;
                  }
                  if (cacheMode == CacheMode::Semantic)
                  {
                      return termsEquivalent(solver, a.formula, b.formula);
                  }
                  if (cacheMode == CacheMode::Canonical)
                  {
                      return solver.simplify(a.formula) == solver.simplify(b.formula);
                  }
                  return false;
              })
{
    d_totalNodes = std::pow(2.0, static_cast<double>(d_projVars.size()));
    d_nodesVisited = 0;
    std::unordered_set<cvc5::Term> all;
    std::function<void(const cvc5::Term&)> collect = [&](const cvc5::Term& t) {
        if (t.getNumChildren() == 0 && t.hasSymbol())
        {
            if (t.toString() != "true" && t.toString() != "false")
            {
                all.insert(t);
            }
            return;
        }
        for (size_t i = 0, n = t.getNumChildren(); i < n; ++i)
        {
            collect(t[i]);
        }
    };
    collect(d_formula);
    std::size_t id = 0;
    for (const auto& v : all)
    {
        d_varToId.emplace(v, id++);
    }
    d_varCount = id;
}

std::uint64_t ProjDDNNF::count()
{
    d_nodesVisited = 0;
    d_cache.clear();
    Profile.resetSearchStats();
    Trace("ddnnf") << "ProjDDNNF count start with " << d_projVars.size()
                   << " variables" << std::endl;
    std::uint64_t res = 0;
    std::vector<std::pair<std::size_t, bool>> assigns;
    if (d_useAssumptions)
    {
        d_solver.push();
        d_solver.assertFormula(d_formula);
        res = countRec(d_formula, d_projVars, 0, assigns);
        d_solver.pop();
    }
    else
    {
        res = countRec(d_formula, d_projVars, 0, assigns);
    }
    Profile.setCacheMemory(d_cache.memoryUsage());
    return res;
}

void ProjDDNNF::collectVars(const cvc5::Term& term,
                            std::unordered_set<std::string>& vars)
{
    if (term.getNumChildren() == 0 && term.hasSymbol())
    {
        std::string name = term.toString();
        if (name != "true" && name != "false")
        {
            vars.insert(name);
        }
        return;
    }
    for (size_t i = 0, n = term.getNumChildren(); i < n; ++i)
    {
        collectVars(term[i], vars);
    }
}

void ProjDDNNF::collectVars(const cvc5::Term& term, sspp::Bitset& vars)
{
    if (term.getNumChildren() == 0 && term.hasSymbol())
    {
        auto it = d_varToId.find(term);
        if (it != d_varToId.end())
        {
            vars.SetTrue(it->second);
        }
        return;
    }
    for (size_t i = 0, n = term.getNumChildren(); i < n; ++i)
    {
        collectVars(term[i], vars);
    }
}

cvc5::Term ProjDDNNF::assignAndRewrite(const cvc5::Term& formula,
                                       const cvc5::Term& var,
                                       bool value,
                                       std::size_t depth)
{
    auto& tm = ttc::getTermBuilder(d_solver);
    cvc5::Term val = tm.mkBoolean(value);
    cvc5::Term substituted = formula.substitute(var, val);
    if (d_propatK > 0 && depth % d_propatK == 0)
    {
        auto start = std::chrono::steady_clock::now();
        cvc5::Term simplified = d_simplifyProp ? d_solver.simplify(substituted)
                                               : simpleRewrite(d_solver, substituted);
        auto end = std::chrono::steady_clock::now();
        Profile.addRewrite(std::chrono::duration<double>(end - start).count());
        return simplified;
    }
    return substituted;
}

void ProjDDNNF::addMonotoneClause()
{
    if (d_falseStack.empty())
    {
        return;
    }
    auto& tm = ttc::getTermBuilder(d_solver);
    cvc5::Term clause = d_falseStack[0];
    for (std::size_t i = 1; i < d_falseStack.size(); ++i)
    {
        clause = tm.mkTerm(cvc5::Kind::OR, {clause, d_falseStack[i]});
    }
    d_solver.assertFormula(clause);
    d_formula = tm.mkTerm(cvc5::Kind::AND, {d_formula, clause});
    d_cache.clear();
    Profile.addMonotoneClause();
}

std::vector<std::size_t> ProjDDNNF::extractConstraintIndices(const cvc5::Term& formula)
{
    std::vector<std::size_t> indices;
    if (formula.getKind() == cvc5::Kind::AND)
    {
        for (size_t i = 0, n = formula.getNumChildren(); i < n; ++i)
        {
            auto s = d_solver.simplify(formula[i]);
            indices.push_back(std::hash<std::string>{}(s.toString()));
        }
    }
    else
    {
        auto s = d_solver.simplify(formula);
        indices.push_back(std::hash<std::string>{}(s.toString()));
    }
    std::sort(indices.begin(), indices.end());
    return indices;
}

std::uint64_t ProjDDNNF::countRec(const cvc5::Term& formula,
                                  const std::vector<cvc5::Term>& projVars,
                                  std::size_t depth,
                                  std::vector<std::pair<std::size_t, bool>>& assigns)
{
    auto& tm = ttc::getTermBuilder(d_solver);
    const bool traceDdnnf = Trace.isEnabled("ddnnf");
    // component cache lookup
    DDNNFCacheKey key;
    if (d_cacheMode == CacheMode::Bool)
    {
        key.assignments = assigns;
        std::sort(key.assignments.begin(), key.assignments.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        key.constraints = extractConstraintIndices(formula);
    }
    else
    {
        cvc5::Term keyFormula =
            (d_cacheMode == CacheMode::Canonical) ? d_solver.simplify(formula) : formula;
        key.formula = keyFormula;
        key.vars = projVars;
    }
    auto startLookup = std::chrono::steady_clock::now();
    Profile.addCacheLookup();
    std::uint64_t cachedVal = 0;
    bool found = d_cache.lookup(key, cachedVal);
    auto endLookup = std::chrono::steady_clock::now();
    Profile.addCacheTime(std::chrono::duration<double>(endLookup - startLookup).count());
    if (found)
    {
        Profile.addCacheHit();
        return cachedVal;
    }

    ++d_nodesVisited;
    double percent = std::min(100.0, (d_nodesVisited / d_totalNodes) * 100.0);
    print_progress(percent, d_cache.size());
    std::unordered_set<std::string> curVars;
    std::size_t totalVarCount = 0;
    if (traceDdnnf)
    {
        collectVars(formula, curVars);
        totalVarCount = curVars.size();
    }
    Trace("ddnnf") << "DDNNF recursion with " << projVars.size()
                   << " projected vars, " << totalVarCount << " total vars"
                   << std::endl;
    Trace("ddnnf") << "Decision level "
                   << (d_projVars.size() - projVars.size()) << std::endl;
    // Attempt decomposition into independent parts
    std::vector<cvc5::Term> conjuncts;
    if (formula.getKind() == cvc5::Kind::AND)
    {
        for (size_t i = 0, n = formula.getNumChildren(); i < n; ++i)
        {
            conjuncts.push_back(formula[i]);
        }
    }
    else
    {
        conjuncts.push_back(formula);
    }

    if (conjuncts.size() > 1)
    {
        size_t n = conjuncts.size();
        if (d_useBitset)
        {
            std::vector<sspp::Bitset> vars(n, sspp::Bitset(d_varCount));
            for (size_t i = 0; i < n; ++i)
            {
                collectVars(conjuncts[i], vars[i]);
            }
            std::vector<size_t> parent(n);
            std::iota(parent.begin(), parent.end(), 0);
            std::function<size_t(size_t)> find = [&](size_t x) {
                if (parent[x] != x) parent[x] = find(parent[x]);
                return parent[x];
            };
            auto unite = [&](size_t a, size_t b) {
                a = find(a); b = find(b);
                if (a != b) parent[b] = a;
            };
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = i + 1; j < n; ++j)
                {
                    if (!(vars[i] & vars[j]).IsEmpty())
                    {
                        unite(i, j);
                    }
                }
            }
            std::unordered_map<size_t, std::vector<size_t>> groups;
            for (size_t i = 0; i < n; ++i)
            {
                groups[find(i)].push_back(i);
            }
            if (groups.size() > 1)
            {
                Trace("decompose") << "decompose " << groups.size() << std::endl;
                Profile.addDecomposition();
                std::uint64_t total = 1;
                for (const auto& kv : groups)
                {
                    const std::vector<size_t>& idxs = kv.second;
                    std::vector<cvc5::Term> groupConj;
                    sspp::Bitset gvars(d_varCount);
                    for (size_t idx : idxs)
                    {
                        groupConj.push_back(conjuncts[idx]);
                        sspp::Bitset tmp = gvars | vars[idx];
                        gvars = std::move(tmp);
                    }
                    cvc5::Term sub = groupConj.size() == 1
                                     ? groupConj[0]
                                     : tm.mkTerm(cvc5::Kind::AND, groupConj);
                    std::vector<cvc5::Term> subProj;
                    for (const auto& pv : projVars)
                    {
                        auto it = d_varToId.find(pv);
                        if (it != d_varToId.end() && gvars.Get(it->second))
                        {
                            subProj.push_back(pv);
                        }
                    }
                    total *= countRec(sub, subProj, depth, assigns);
                }
                d_cache.insert(key, total);
                return total;
            }
        }
        else
        {
            std::vector<std::unordered_set<std::string>> vars(n);
            for (size_t i = 0; i < n; ++i)
            {
                collectVars(conjuncts[i], vars[i]);
            }
            std::vector<size_t> parent(n);
            std::iota(parent.begin(), parent.end(), 0);
            std::function<size_t(size_t)> find = [&](size_t x) {
                if (parent[x] != x) parent[x] = find(parent[x]);
                return parent[x];
            };
            auto unite = [&](size_t a, size_t b) {
                a = find(a); b = find(b);
                if (a != b) parent[b] = a;
            };
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = i + 1; j < n; ++j)
                {
                    bool share = false;
                    for (const auto& v : vars[i])
                    {
                        if (vars[j].count(v))
                        {
                            share = true; break;
                        }
                    }
                    if (share)
                    {
                        unite(i, j);
                    }
                }
            }
            std::unordered_map<size_t, std::vector<size_t>> groups;
            for (size_t i = 0; i < n; ++i)
            {
                groups[find(i)].push_back(i);
            }
            if (groups.size() > 1)
            {
                Trace("decompose") << "decompose " << groups.size() << std::endl;
                Profile.addDecomposition();
                std::uint64_t total = 1;
                for (const auto& kv : groups)
                {
                    const std::vector<size_t>& idxs = kv.second;
                    std::vector<cvc5::Term> groupConj;
                    std::unordered_set<std::string> gvars;
                    for (size_t idx : idxs)
                    {
                        groupConj.push_back(conjuncts[idx]);
                        gvars.insert(vars[idx].begin(), vars[idx].end());
                    }
                    cvc5::Term sub = groupConj.size() == 1
                                     ? groupConj[0]
                                     : tm.mkTerm(cvc5::Kind::AND, groupConj);
                    std::vector<cvc5::Term> subProj;
                    for (const auto& pv : projVars)
                    {
                        std::string name = pv.toString();
                        if (gvars.count(name))
                        {
                            subProj.push_back(pv);
                        }
                    }
                    total *= countRec(sub, subProj, depth, assigns);
                }
                d_cache.insert(key, total);
                return total;
            }
        }
    }

    // Branch on projection variables when no decomposition is possible
    if (projVars.empty())
    {
        std::uint64_t res = 0;
        if (d_useAssumptions)
        {
            Profile.addCheckSatCall();
            cvc5::Result r = d_solver.checkSatAssuming(d_assumptions);
            if (r.isUnsat()) Profile.addUnsatCheckSatCall();
            res = r.isSat() ? 1 : 0;
        }
        else
        {
            d_solver.push();
            d_solver.assertFormula(formula);
            Profile.addCheckSatCall();
            cvc5::Result r = d_solver.checkSat();
            if (r.isUnsat()) Profile.addUnsatCheckSatCall();
            d_solver.pop();
            res = r.isSat() ? 1 : 0;
        }
        if (d_useMono && res == 0)
        {
            addMonotoneClause();
        }
        d_cache.insert(key, res);
        return res;
    }

    const cvc5::Term& var = projVars[0];
    std::vector<cvc5::Term> rest(projVars.begin() + 1, projVars.end());
    std::uint64_t total = 0;
    std::size_t newDepth = depth + 1;

    std::array<bool, 2> order;
    if (d_useMono)
    {
        bool def = d_monoTrue;
        order = {def, !def};
    }
    else
    {
        order = {true, false};
    }

    auto itId = d_varToId.find(var);
    std::size_t vid = (itId != d_varToId.end()) ? itId->second : 0;
    for (bool val : order)
    {
        Profile.addDecision();
        Trace("ddnnf") << "Decision " << var << " = " << (val ? "true" : "false") << std::endl;
        cvc5::Term branch = assignAndRewrite(formula, var, val, newDepth);
        Trace("ddnnf") << "Residual formula: " << branch << std::endl;
        bool pushFalse = d_useMono && !val;
        if (pushFalse)
        {
            d_falseStack.push_back(var);
        }

        if (d_useAssumptions)
        {
            cvc5::Term lit = val ? var : tm.mkTerm(cvc5::Kind::NOT, {var});
            d_assumptions.push_back(lit);
            Profile.addCheckSatCall();
            cvc5::Result r = d_solver.checkSatAssuming(d_assumptions);
            if (r.isUnsat())
            {
                Profile.addUnsatCheckSatCall();
            }
            if (!r.isUnsat())
            {
                assigns.emplace_back(vid, val);
                total += countRec(branch, rest, newDepth, assigns);
                assigns.pop_back();
            }
            else if (d_useMono)
            {
                addMonotoneClause();
            }
            d_assumptions.pop_back();
        }
        else
        {
            d_solver.push();
            d_solver.assertFormula(branch);
            Profile.addCheckSatCall();
            cvc5::Result r = d_solver.checkSat();
            if (r.isUnsat())
            {
                Profile.addUnsatCheckSatCall();
            }
            if (!r.isUnsat())
            {
                assigns.emplace_back(vid, val);
                total += countRec(branch, rest, newDepth, assigns);
                assigns.pop_back();
            }
            else if (d_useMono)
            {
                d_solver.pop();
                addMonotoneClause();
                if (pushFalse)
                {
                    d_falseStack.pop_back();
                }
                continue;
            }
            d_solver.pop();
        }

        if (pushFalse)
        {
            d_falseStack.pop_back();
        }
    }

    d_cache.insert(key, total);
    return total;
}
