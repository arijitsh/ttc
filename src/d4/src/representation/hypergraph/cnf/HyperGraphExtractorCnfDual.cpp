/*
 * d4
 * Copyright (C) 2020  Univ. Artois & CNRS
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#include "HyperGraphExtractorCnfDual.hpp"

#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "src/formulaManager/cnf/CnfManager.hpp"
#include "../../../../../d4_smt.hpp"

namespace d4 {

/**
 * @brief HyperGraphExtractorCnfDual::~HyperGraphExtractorCnfDual
 * implementation.
 */
HyperGraphExtractorCnfDual::~HyperGraphExtractorCnfDual() {
  if (m_data) delete[] m_data;
}  // destructor

/**
 * @brief HyperGraphExtractorCnfDual::constructHyperGraph implementation.
 */
namespace {
// Collect all SMT-level variables occurring in a term.
void collectVars(const cvc5::Term &t, std::unordered_set<std::string> &vars) {
  if (t.getNumChildren() == 0 && t.hasSymbol()) {
    std::string name = t.toString();
    if (name != "true" && name != "false") vars.insert(name);
    return;
  }
  for (size_t i = 0, n = t.getNumChildren(); i < n; ++i)
    collectVars(t[i], vars);
}

int find(std::vector<int> &parent, int x) {
  if (parent[x] != x) parent[x] = find(parent, parent[x]);
  return parent[x];
}

void unite(std::vector<int> &parent, int a, int b) {
  a = find(parent, a);
  b = find(parent, b);
  if (a != b) parent[b] = a;
}
}  // namespace

InfoHyperGraph HyperGraphExtractorCnfDual::constructHyperGraph(
    FormulaManager &om, std::vector<Var> &component, HyperGraph &hypergraph) {
  CnfManager &tmp = static_cast<CnfManager &>(om);

  // Union-find over SMT variables occurring in CNF variables
  std::unordered_map<std::string, int> nameToIdx;
  std::vector<int> parent;
  std::unordered_map<Var, std::vector<std::string>> varNames;
  for (Var v : component) {
    if (v >= (int)d4_smt_idxToTerm.size()) continue;
    const cvc5::Term &t = d4_smt_idxToTerm[v];
    std::unordered_set<std::string> vars;
    if (!t.isNull()) collectVars(t, vars);
    std::vector<int> ids;
    ids.reserve(vars.size());
    for (const auto &name : vars) {
      auto [it, inserted] = nameToIdx.emplace(name, parent.size());
      if (inserted) parent.push_back(it->second);
      ids.push_back(it->second);
    }
    if (!ids.empty()) {
      for (size_t i = 1; i < ids.size(); ++i) unite(parent, ids[0], ids[i]);
    }
    varNames.emplace(v, std::vector<std::string>(vars.begin(), vars.end()));
  }

  // Group CNF variables by the SMT variables they mention
  std::unordered_map<int, Var> rootToRep;
  m_groups.clear();
  for (Var v : component) {
    const auto &names = varNames[v];
    if (names.empty()) {
      m_groups[v].push_back(v);
      continue;
    }
    int idx = nameToIdx[names[0]];
    int root = find(parent, idx);
    Var rep;
    auto it = rootToRep.find(root);
    if (it == rootToRep.end()) {
      rep = v;
      rootToRep[root] = v;
      m_groups[rep];
    } else {
      rep = it->second;
    }
    m_groups[rep].push_back(v);
  }

  // allocate memory
  unsigned pos = 0;
  m_data = new char[component.size() * sizeof(HyperEdge) +
                    sizeof(unsigned) * tmp.getSumSizeClauses()];

  // create the graph with grouped variables
  for (auto &p : m_groups) {
    Var rep = p.first;
    std::vector<unsigned> edgeData;
    for (Var v : p.second) {
      for (auto &l : {Lit::makeLitFalse(v), Lit::makeLitTrue(v)}) {
        IteratorIdxClause listIdx = tmp.getVecIdxClause(l);
        for (int *ptr = listIdx.start; ptr != listIdx.end; ++ptr)
          edgeData.push_back(*ptr);
      }
    }
    hypergraph.addEdge(new (&m_data[pos])
                           HyperEdge((unsigned)rep, edgeData.size(),
                                     edgeData.data()));
    pos += sizeof(HyperEdge) + edgeData.size() * sizeof(unsigned);
  }

  return {dynamic_cast<CnfManager &>(om).getNbVariable(),
          dynamic_cast<CnfManager &>(om).getNbClause(),
          dynamic_cast<CnfManager &>(om).getSumSizeClauses()};
}  // constructHyperGraph

/**
 * @brief HyperGraphExtractorCnfDual::split implementation.
 */
void HyperGraphExtractorCnfDual::split(HyperGraph &graph,
                                       std::vector<int> &partition,
                                       std::vector<Var> &cut,
                                       HyperGraph &firstGraph,
                                       HyperGraph &secondGraph) {
  auto pushGroup = [&](Var rep) {
    auto it = m_groups.find(rep);
    if (it != m_groups.end())
      cut.insert(cut.end(), it->second.begin(), it->second.end());
    else
      cut.push_back(rep);
  };

  if (graph.getNbEdges() < 10) {
    for (unsigned i = 0; i < graph.getNbEdges(); i++) pushGroup(graph[i].getId());
    return;
  }

  for (unsigned i = 0; i < graph.getNbEdges(); i++) {
    HyperEdge &e = graph[i];
    if (!e.getSize()) continue;

    int part = partition[e[0]];
    bool clash = false;
    for (unsigned j = 1; !clash && j < e.getSize(); j++)
      clash = part != partition[e[j]];

    if (clash)
      pushGroup(e.getId());
    else if (part == 0)
      firstGraph.addEdge(graph.getEdge(i));
    else
      secondGraph.addEdge(graph.getEdge(i));
  }

  if (!firstGraph.getNbEdges()) {
    secondGraph.setNbEdges(0);
    for (unsigned i = 0; i < graph.getNbEdges(); i++) pushGroup(graph[i].getId());
  }

  if (!secondGraph.getNbEdges()) {
    firstGraph.setNbEdges(0);
    for (unsigned i = 0; i < graph.getNbEdges(); i++) pushGroup(graph[i].getId());
  }
}  // split

}  // namespace d4
