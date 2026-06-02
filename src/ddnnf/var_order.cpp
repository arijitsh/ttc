#include "var_order.hpp"

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <limits>
#include "logger.hpp"
#include "profiler.hpp"

// FlowCutter wrapper
#include "iflowcutter.hpp"
// Graph structure and contraction
#include "TreeDecomposition.hpp"

namespace {

// Collect all variable names occurring in a term
void collectVars(const cvc5::Term& term, std::unordered_set<std::string>& vars)
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

} // namespace

// Compute the treewidth of the variable incidence graph given an elimination order
int computeTreewidth(const ArrayIDIDFunc& order,
                     int nodeCount,
                     const std::vector<int>& tails,
                     const std::vector<int>& heads)
{
    std::vector<std::unordered_set<int>> adj(nodeCount);
    for (size_t i = 0; i < tails.size(); ++i)
    {
        int u = tails[i];
        int v = heads[i];
        if (u == v) continue;
        adj[u].insert(v);
        adj[v].insert(u);
    }

    int tw = 0;
    for (int idx = 0; idx < nodeCount; ++idx)
    {
        int x = order[idx];
        std::vector<int> neigh(adj[x].begin(), adj[x].end());
        if (static_cast<int>(neigh.size()) > tw)
        {
            tw = neigh.size();
        }
        for (int u : neigh)
        {
            for (int v : neigh)
            {
                if (u != v)
                {
                    adj[u].insert(v);
                    adj[v].insert(u);
                }
            }
        }
        for (int u : neigh)
        {
            adj[u].erase(x);
        }
        adj[x].clear();
    }
    return tw;
}

std::vector<cvc5::Term> computeProjVarOrder(const cvc5::Term& formula,
                                            const std::vector<cvc5::Term>& projVars,
                                            cvc5::Solver& solver,
                                            bool printTD,
                                            bool contract,
                                            bool netrel)
{
    (void)solver; // unused in the heuristic for now

    if (projVars.size() <= 1)
    {
        return projVars;
    }

    // Break the formula into top-level assertions/conjuncts
    std::vector<cvc5::Term> asserts;
    if (formula.getKind() == cvc5::Kind::AND)
    {
        for (size_t i = 0, n = formula.getNumChildren(); i < n; ++i)
        {
            asserts.push_back(formula[i]);
        }
    }
    else
    {
        asserts.push_back(formula);
    }

    // First collect all variable names in the formula to assign ids
    std::unordered_set<std::string> allVars;
    collectVars(formula, allVars);
    std::unordered_map<std::string, int> varToId;
    int nextId = 0;
    for (const auto& v : allVars)
    {
        varToId.emplace(v, nextId++);
    }

    // Build primal graph on all variables
    TWD::Graph primal(nextId);
    for (const auto& as : asserts)
    {
        std::unordered_set<std::string> vars;
        collectVars(as, vars);
        std::vector<int> ids;
        ids.reserve(vars.size());
        for (const auto& v : vars)
        {
            ids.push_back(varToId[v]);
        }
        for (size_t i = 0; i < ids.size(); ++i)
        {
            for (size_t j = i + 1; j < ids.size(); ++j)
            {
                primal.addEdge(ids[i], ids[j]);
            }
        }
    }

    // Determine ids of projection variables
    std::vector<int> projIds;
    projIds.reserve(projVars.size());
    std::unordered_map<int, cvc5::Term> idToProjTerm;
    std::vector<cvc5::Term> missing;
    std::vector<int> newToOld;
    for (const auto& v : projVars)
    {
        std::string name = v.toString();
        auto it = varToId.find(name);
        if (it != varToId.end())
        {
            projIds.push_back(it->second);
            idToProjTerm.emplace(it->second, v);
        }
        else
        {
            missing.push_back(v);
        }
    }

    std::size_t merged = 0;
    if (netrel)
    {
        std::vector<int> isProj(nextId, 0);
        for (int id : projIds) isProj[id] = 1;
        std::vector<int> newProjIds;
        std::unordered_map<int, cvc5::Term> newIdToProjTerm;
        for (int pid : projIds)
        {
            const auto& neigh = primal.Neighbors(pid);
            if (neigh.size() == 1)
            {
                int neighId = neigh[0];
                if (!isProj[neighId])
                {
                    primal.contract(pid, std::numeric_limits<int>::max());
                    newProjIds.push_back(neighId);
                    newIdToProjTerm.emplace(neighId, idToProjTerm[pid]);
                    isProj[neighId] = 1;
                    ++merged;
                    continue;
                }
            }
            newProjIds.push_back(pid);
            newIdToProjTerm.emplace(pid, idToProjTerm[pid]);
        }
        projIds.swap(newProjIds);
        idToProjTerm.swap(newIdToProjTerm);
    }

    if (projIds.size() <= 1)
    {
        return projVars;
    }

    // Record counts before contraction
    std::size_t preNodes = nextId - merged;
    std::size_t preEdges = primal.numEdges();

    TWD::IFlowCutter fc;
    std::vector<int> tails;
    std::vector<int> heads;
    ArrayIDIDFunc order;
    int graphNodes = 0;
    double decompTime = 0.0;
    int tw = 0;

    if (contract)
    {
        // Contract variables not in the projection set
        std::vector<int> isProj(nextId, 0);
        for (int id : projIds) isProj[id] = 1;
        for (int i = 0; i < nextId; ++i)
        {
            if (!isProj[i])
            {
                primal.contract(i, std::numeric_limits<int>::max());
            }
        }

        // Build graph containing only projection variables
        TWD::Graph projGraph(static_cast<int>(projIds.size()));
        std::unordered_map<int, int> oldToNew;
        newToOld.resize(projIds.size());
        for (size_t i = 0; i < projIds.size(); ++i)
        {
            oldToNew.emplace(projIds[i], static_cast<int>(i));
            newToOld[i] = projIds[i];
        }
        for (size_t i = 0; i < projIds.size(); ++i)
        {
            int u_orig = projIds[i];
            for (int v_orig : primal.Neighbors(u_orig))
            {
                if (isProj[v_orig] && u_orig < v_orig)
                {
                    int u_new = static_cast<int>(i);
                    int v_new = oldToNew[v_orig];
                    projGraph.addEdge(u_new, v_new);
                }
            }
        }

        const auto& adj = projGraph.get_adj_list();
        for (int u = 0; u < projGraph.numNodes(); ++u)
        {
            for (int v : adj[u])
            {
                tails.push_back(u);
                heads.push_back(v);
            }
        }

        auto start = std::chrono::steady_clock::now();
        fc.importGraph(projGraph.numNodes(), tails, heads);
        order = fc.constructTD();
        auto end = std::chrono::steady_clock::now();
        decompTime = std::chrono::duration<double>(end - start).count();

        graphNodes = projGraph.numNodes();
        tw = computeTreewidth(order, projGraph.numNodes(), tails, heads);

        if (printTD)
        {
            std::size_t edgeCount = tails.size() / 2;
            double density = projGraph.numNodes() > 1 ? static_cast<double>(edgeCount) / (projGraph.numNodes() * (projGraph.numNodes() - 1) / 2.0) : 0.0;
            double edgeVar = projVars.empty() ? 0.0 : static_cast<double>(edgeCount) / projVars.size();
            print_section("tree decomposition");
            if (netrel)
            {
                std::cout << "c after netrel nodes: " << preNodes << " edges: " << preEdges << std::endl;
            }
            else
            {
                std::cout << "c before contraction nodes: " << preNodes << " edges: " << preEdges << std::endl;
            }
            std::cout << "c nodes: " << projGraph.numNodes() << " nvars: " << projVars.size()
                      << " edges: " << edgeCount << std::endl;
            std::cout << "c Primal graph   nodes: " << projGraph.numNodes() << " edges: " << edgeCount
                      << " density: " << std::fixed << std::setprecision(3) << density
                      << " edge/var: " << std::fixed << std::setprecision(3) << edgeVar << std::endl;
            std::cout << "c min degree heuristic" << std::endl;
            std::cout << "c #bags " << (projGraph.numNodes() > 0 ? projGraph.numNodes() - 1 : 0)
                      << ", tw " << tw << ", elapsed 0.000 s" << std::endl;
            std::cout << "c min shortcut heuristic" << std::endl;
            std::cout << "c decompose time: " << std::fixed << std::setprecision(3)
                      << decompTime << std::endl;
            std::cout << "c" << std::endl;
        }

    }
    else
    {
        const auto& adj = primal.get_adj_list();
        for (int u = 0; u < primal.numNodes(); ++u)
        {
            for (int v : adj[u])
            {
                tails.push_back(u);
                heads.push_back(v);
            }
        }

        auto start = std::chrono::steady_clock::now();
        fc.importGraph(primal.numNodes(), tails, heads);
        order = fc.constructTD();
        auto end = std::chrono::steady_clock::now();
        decompTime = std::chrono::duration<double>(end - start).count();

        graphNodes = primal.numNodes();
        tw = computeTreewidth(order, primal.numNodes(), tails, heads);

        if (printTD)
        {
            std::size_t edgeCount = tails.size() / 2;
            double density = primal.numNodes() > 1 ? static_cast<double>(edgeCount) / (primal.numNodes() * (primal.numNodes() - 1) / 2.0) : 0.0;
            double edgeVar = projVars.empty() ? 0.0 : static_cast<double>(edgeCount) / projVars.size();
            print_section("tree decomposition");
            if (netrel)
            {
                std::cout << "c after netrel nodes: " << preNodes << " edges: " << preEdges << std::endl;
            }
            std::cout << "c nodes: " << primal.numNodes() << " nvars: " << projVars.size()
                      << " edges: " << edgeCount << std::endl;
            std::cout << "c Primal graph   nodes: " << primal.numNodes() << " edges: " << edgeCount
                      << " density: " << std::fixed << std::setprecision(3) << density
                      << " edge/var: " << std::fixed << std::setprecision(3) << edgeVar << std::endl;
            std::cout << "c min degree heuristic" << std::endl;
            std::cout << "c #bags " << (primal.numNodes() > 0 ? primal.numNodes() - 1 : 0)
                      << ", tw " << tw << ", elapsed 0.000 s" << std::endl;
            std::cout << "c min shortcut heuristic" << std::endl;
            std::cout << "c decompose time: " << std::fixed << std::setprecision(3)
                      << decompTime << std::endl;
            std::cout << "c" << std::endl;
        }
    }

    const bool traceVarOrder = Trace.isEnabled("var-order");
    if (traceVarOrder)
    {
        Trace("var-order") << "Treewidth: " << tw << std::endl;
    }

    // Extract ordered projection variables
    std::vector<cvc5::Term> ordered;
    ordered.reserve(projIds.size() + missing.size());
    for (int i = 0; i < graphNodes; ++i)
    {
        int node = order[i];
        if (contract)
        {
            node = newToOld[node];
        }
        auto it = idToProjTerm.find(node);
        if (it != idToProjTerm.end())
        {
            ordered.push_back(it->second);
        }
    }
    ordered.insert(ordered.end(), missing.begin(), missing.end());

    if (traceVarOrder)
    {
        auto stream = Trace("var-order");
        stream << "Variable order: ";
        for (const auto& v : ordered)
        {
            stream << v.toString() << ' ';
        }
        stream << std::endl;
    }
    Profile.addTreeDecomp(decompTime);
    return ordered;
}
