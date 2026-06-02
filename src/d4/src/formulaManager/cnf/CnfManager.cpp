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

#include "CnfManager.hpp"

#include <algorithm>  // std::sort
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "CnfManagerDyn.hpp"
#include "ResidualSimplifier.hpp"
#include "../../../../d4_smt.hpp"
#include "../../../../ddnnf/logger.hpp"
#include "../../../../ddnnf/profiler.hpp"
#include "src/methods/nnf/Node.hpp"
#include "src/problem/ProblemTypes.hpp"

namespace d4 {

namespace {
// Collect SMT-level variable names occurring in a term (leaf symbols only).
inline void collectSmtVars(const cvc5::Term &t,
                           std::unordered_set<std::string> &vars) {
  if (t.getNumChildren() == 0 && t.hasSymbol()) {
    const std::string name = t.getSymbol();
    if (name != "true" && name != "false") vars.insert(name);
    return;
  }
  const size_t n = t.getNumChildren();
  for (size_t i = 0; i < n; ++i) collectSmtVars(t[i], vars);
}

// Cache, per CNF variable index, the set of SMT variable names that occur
// in its associated term. This avoids re-traversing the SMT term at every
// recursive call while computing connected components.
static std::vector<std::vector<std::string>> smtVarNamesCache;

inline void ensureSmtNamesCache() {
  const size_t n = d4_smt_idxToTerm.size();
  if (smtVarNamesCache.size() == n) return;
  smtVarNamesCache.clear();
  smtVarNamesCache.resize(n);
  for (size_t i = 0; i < n; ++i) {
    const cvc5::Term &t = d4_smt_idxToTerm[i];
    if (t.isNull()) continue;
    std::unordered_set<std::string> names;
    collectSmtVars(t, names);
    smtVarNamesCache[i].assign(names.begin(), names.end());
  }
}
}  // namespace

/**
   Constructor.
*/
CnfManager::CnfManager(ProblemManager &p) : FormulaManager(p.getNbVar()) {
  // get the clauses.
  try {
    CnfMatrix &pcnf = dynamic_cast<CnfMatrix &>(p);
    m_clauses = pcnf.getClauses();
  } catch (std::bad_cast &bc) {
    std::cerr << "c bad_cast caught: " << bc.what() << '\n';
    std::cerr << "c A CNF formula was expeted\n";
    assert(0);
  }

  // store the not binary clauses.
  m_maxSizeClause = 0;
  unsigned count = 0;
  std::vector<std::vector<int>> occurrence((m_nbVar + 1) << 1);
  for (unsigned i = 0; i < m_clauses.size(); i++) {
    if (m_clauses[i].size() > 2) m_clausesNotBin.push_back(i);

    for (auto &l : m_clauses[i]) occurrence[l.intern()].push_back(i);
    count += m_clauses[i].size();

    if (m_clauses[i].size() > m_maxSizeClause)
      m_maxSizeClause = m_clauses[i].size();
  }

  // reserve the memory to store the occurrence lists.
  m_occurrence.resize((m_nbVar + 1) << 1, {NULL, NULL, 0, 0});
  m_dataOccurrenceMemory = new int[count];

  // construct the occurrence list.
  int *ptr = m_dataOccurrenceMemory;
  for (unsigned i = 0; i < occurrence.size(); i++) {
    std::vector<int> &occList = occurrence[i];

    unsigned posNotBin = occList.size() - 1;
    for (auto const &idx : occList) {
      if (m_clauses[idx].size() == 2)
        ptr[m_occurrence[i].nbBin++] = idx;
      else
        ptr[posNotBin--] = idx;
    }

    m_occurrence[i].bin = ptr;
    m_occurrence[i].notBin = &ptr[posNotBin + 1];
    m_occurrence[i].nbNotBin = occList.size() - m_occurrence[i].nbBin;
    ptr = &ptr[occList.size()];
  }

  // variables:
  m_inCurrentComponent.resize(m_nbVar + 1, false);
  m_idxComponent.resize(m_nbVar + 1, 0);

  // clauses:
  unsigned nbClause = m_clauses.size();
  m_mustUnMark.reserve(nbClause);
  m_markView.assign(nbClause, 0);

  m_infoClauses.resize(nbClause);

  // set the info about xorLitBin.
  for (unsigned i = 0; i < m_clauses.size(); i++) {
    if (m_clauses[i].size() == 2)
      m_infoClauses[i].xorLitBin =
          m_clauses[i][0].intern() ^ m_clauses[i][1].intern();
  }

  m_infoCluster.resize(p.getNbVar() + nbClause + 1, {0, 0, -1});
  m_activeVariables = new Var[p.getNbVar() + 1];

  m_occInitSizeNotBin.resize((p.getNbVar() + 1) << 1, 0);
  assert(m_occInitSizeNotBin.size() == m_occurrence.size());
  for (unsigned i = 0; i < m_occurrence.size(); i++)
    m_occInitSizeNotBin[i] = m_occurrence[i].nbNotBin;
}  // construtor

/**
 * @brief Destroy the Spec Manager Cnf:: Spec Manager Cnf object
 *
 */
CnfManager::~CnfManager() { delete[] m_dataOccurrenceMemory; }  // destructor

/**
   Look all the formula in order to compute the connected component
   of the formula (union find algorithm).

   @param[out] varCo, the different connected components found
   @param[in] setOfVar, the current set of variables
   @param[out] freeVar, the set of variables that are present in setOfVar but
   not in the problem anymore

   \return the number of component found
*/
int CnfManager::computeConnectedComponent(std::vector<std::vector<Var>> &varCo,
                                          std::vector<Var> &setOfVar,
                                          std::vector<Var> &freeVar) {
  const bool traceDecompose = Trace.isEnabled("decompose");
  Var *lastVar = m_activeVariables, *currentVar = m_activeVariables;
  for (auto v : setOfVar) {
    assert(v < m_infoCluster.size());
    m_infoCluster[v].parent = v;
    m_infoCluster[v].size = 1;
    if (m_currentValue[v] == l_Undef) *lastVar++ = v;
  }

  while (currentVar != lastVar) {
    Var v = *currentVar++;
    assert(m_currentValue[v] == l_Undef);

    // visit the index clauses
    Var rootV = findRoot(v);
    Lit l = Lit::makeLit(v, false);

    for (unsigned i = 0; i < 2; i++) {  // both literals.
      IteratorIdxClause listIndex = getVecIdxClause(l);
      for (int *ptr = listIndex.start; ptr != listIndex.end; ptr++) {
        int idx = *ptr;
        // Guard against out-of-range clause indices that could otherwise lead
        // to undefined behaviour when accessing m_markView.  Such indices may
        // appear if the occurrence list becomes inconsistent.  In that case we
        // simply skip the clause since it cannot contribute to any component.
        if (idx < 0 || static_cast<size_t>(idx) >= m_markView.size()) {
          continue;
        }

        if (m_infoClauses[idx].isSat) continue;

        const Var clauseNode = idx + m_nbVar + 1;

        if (!m_markView[idx]) {
          m_markView[idx] = 1;
          m_infoCluster[clauseNode].parent = rootV;
          m_infoCluster[rootV].size++;
          m_mustUnMark.push_back(idx);
          continue;
        }

        Var rootW = findRoot(m_infoCluster[clauseNode].parent);
        rootV = findRoot(rootV);
        if (rootV == rootW) {
          m_infoCluster[clauseNode].parent = rootV;
          continue;
        }

        if (m_infoCluster[rootV].size < m_infoCluster[rootW].size) {
          m_infoCluster[rootW].size += m_infoCluster[rootV].size;
          m_infoCluster[rootV].parent = rootW;
          rootV = rootW;
        } else {
          m_infoCluster[rootV].size += m_infoCluster[rootW].size;
          m_infoCluster[rootW].parent = rootV;
        }

        m_infoCluster[clauseNode].parent = rootV;
      }

      l = ~l;
    }
  }

  int cnfComponentCount = 0;
  for (Var *it = m_activeVariables; it != lastVar; ++it) {
    Var v = *it;
    if (m_currentValue[v] != l_Undef) continue;
    Var rootV = findRoot(v);
    if (rootV == v && m_infoCluster[v].size > 1) ++cnfComponentCount;
  }
  std::unordered_map<Var, std::vector<Var>> componentTrace;
  auto describeComponent = [&](Var root) -> std::string {
    if (!traceDecompose) return "{}";
    root = findRoot(root);
    auto it = componentTrace.find(root);
    if (it == componentTrace.end() || it->second.empty()) return "{}";
    std::ostringstream ss;
    ss << '{';
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (i) ss << ", ";
      ss << it->second[i];
    }
    ss << '}';
    return ss.str();
  };
  auto mergeComponentData = [&](Var keepRoot, Var otherRoot) {
    if (!traceDecompose) return;
    keepRoot = findRoot(keepRoot);
    otherRoot = findRoot(otherRoot);
    if (keepRoot == otherRoot) return;
    auto itKeep = componentTrace.find(keepRoot);
    auto itOther = componentTrace.find(otherRoot);
    if (itKeep == componentTrace.end()) {
      if (itOther == componentTrace.end()) return;
      componentTrace.emplace(keepRoot, itOther->second);
      componentTrace.erase(itOther);
      return;
    }
    if (itOther == componentTrace.end()) return;
    auto &dest = itKeep->second;
    dest.insert(dest.end(), itOther->second.begin(), itOther->second.end());
    componentTrace.erase(itOther);
  };
  if (traceDecompose) {
    // Group variables by their root component
    std::unordered_map<Var, std::vector<Var>> componentSets;
    for (Var *it = m_activeVariables; it != lastVar; ++it) {
      Var v = *it;
      if (m_currentValue[v] != l_Undef) continue;
      Var rootV = findRoot(v);
      componentSets[rootV].push_back(v);
      componentTrace[rootV].push_back(v);
    }

    // Print components with size > 1
    Trace("decompose") << "c Component variable sets:\n";
    int componentNum = 0;
    for (const auto& component : componentSets) {
      if (component.second.size() > 1) {
        std::ostringstream ss;
        ss << "c Component " << componentNum++ << ": {";
        for (size_t i = 0; i < component.second.size(); i++) {
          ss << component.second[i];
          if (i < component.second.size() - 1) ss << ", ";
        }
        ss << "}";
        Trace("decompose") << ss.str() << "\n";
      }
    }
  }
  if (traceDecompose) {
    Trace("decompose") << "c CNF connected components found: "
                        << cnfComponentCount << "\n";
  }

  // Merge components that share SMT-level variables.
  std::unordered_map<std::string, Var> nameToRoot;
  const std::vector<cvc5::Term> *residual = nullptr;
  const std::vector<std::vector<std::string>> *residualVarLists = nullptr;
  const std::unordered_set<std::string> *residualVarNames = nullptr;
  bool considerSmtDependencies = !d4_full_decompose;

  m_residualAssertions.clear();
  m_residualAssertionVars.clear();
  m_residualVarNames.clear();

  if (considerSmtDependencies && d4_smt_solver) {
    ensureSmtNamesCache();

    auto assertions = d4_smt_solver->getAssertions();
    if (m_residualAssertions.capacity() < assertions.size())
      m_residualAssertions.reserve(assertions.size());

    for (const auto &a : assertions) {
      cvc5::Term r = a;

      if (r.isBooleanValue() && r.getBooleanValue()) continue;

      m_residualAssertions.push_back(r);
    }

    if (d4_use_residual_simplifier && !m_residualAssertions.empty()) {
      ResidualSimplifier simplifier(*d4_smt_solver);
      SimplifyResult simplifyResult = simplifier.simplify(m_residualAssertions);
      m_residualAssertions.swap(simplifyResult.assertions);
    }

    m_residualAssertionVars.clear();
    m_residualVarNames.clear();
    if (!m_residualAssertions.empty()) {
      if (m_residualAssertionVars.capacity() < m_residualAssertions.size())
        m_residualAssertionVars.reserve(m_residualAssertions.size());

      for (const auto &r : m_residualAssertions) {
        m_residualAssertionVars.emplace_back();
        auto &names = m_residualAssertionVars.back();
        if (r.isBooleanValue()) continue;

        m_residualNameScratch.clear();
        collectSmtVars(r, m_residualNameScratch);

        names.assign(m_residualNameScratch.begin(), m_residualNameScratch.end());
        m_residualVarNames.insert(names.begin(), names.end());
      }
    }

    residual = &m_residualAssertions;
    residualVarLists = &m_residualAssertionVars;
    residualVarNames = &m_residualVarNames;

    if (traceDecompose) {
      Trace("decompose") << "c Residual SMT formula clauses: "
                          << m_residualAssertions.size() << "\n";
    }
  }

  if (considerSmtDependencies && !d4_smt_idxToTerm.empty()) {
    ensureSmtNamesCache();
    if (residualVarNames && !residualVarNames->empty()) {
      if (traceDecompose) {
        Trace("decompose")
            << "c Merging components connected via residual SMT variable names\n";
      }

      // Print which CNF variables share SMT variables
      if (traceDecompose) {
        std::unordered_map<std::string, std::vector<Var>> smtVarToCnfVars;
        for (Var *it = m_activeVariables; it != lastVar; ++it) {
          Var v = *it;
          if (v >= (int)d4_smt_idxToTerm.size() ||
              m_currentValue[v] != l_Undef)
            continue;
          const cvc5::Term &t = d4_smt_idxToTerm[v];
          if (t.isNull()) continue;
          const auto &vars = smtVarNamesCache[v];
          for (const auto &name : vars) {
            if (residualVarNames->find(name) == residualVarNames->end())
              continue;
            smtVarToCnfVars[name].push_back(v);
          }
        }

        // Print SMT variables that connect multiple CNF variables
        Trace("decompose") << "c SMT variables connecting CNF variables:\n";
        for (const auto &entry : smtVarToCnfVars) {
          if (entry.second.size() > 1) {
            std::ostringstream ss;
            ss << "c SMT var '" << entry.first << "' connects CNF vars: {";
            for (size_t i = 0; i < entry.second.size(); i++) {
              ss << entry.second[i];
              if (i < entry.second.size() - 1) ss << ", ";
            }
            ss << "}";
            Trace("decompose") << ss.str() << "\n";
          }
        }
      }

      for (Var *it = m_activeVariables; it != lastVar; ++it) {
        Var v = *it;
        if (v >= (int)d4_smt_idxToTerm.size()) continue;
        const cvc5::Term &t = d4_smt_idxToTerm[v];
        if (t.isNull()) continue;
        const auto &vars = smtVarNamesCache[v];
        for (const auto &name : vars) {
          if (residualVarNames->find(name) == residualVarNames->end()) continue;

          Var rootV = findRoot(v);
          auto itName = nameToRoot.find(name);
          if (itName == nameToRoot.end()) {
            nameToRoot.emplace(name, rootV);
            continue;
          }
          Var rootW = findRoot(itName->second);
          if (rootV == rootW) continue;
          if (traceDecompose) {
            auto itCompV = componentTrace.find(rootV);
            auto itCompW = componentTrace.find(rootW);
            if (itCompV != componentTrace.end() && !itCompV->second.empty() &&
                itCompW != componentTrace.end() && !itCompW->second.empty()) {
              Trace("decompose") << "SMT variable '" << name
                                 << "' merges CNF components "
                                 << describeComponent(rootV) << " and "
                                 << describeComponent(rootW) << std::endl;
            }
          }
          if (m_infoCluster[rootV].size < m_infoCluster[rootW].size) {
            if (traceDecompose) mergeComponentData(rootW, rootV);
            m_infoCluster[rootW].size += m_infoCluster[rootV].size;
            m_infoCluster[rootV].parent = rootW;
            rootV = rootW;
          } else {
            if (traceDecompose) mergeComponentData(rootV, rootW);
            m_infoCluster[rootV].size += m_infoCluster[rootW].size;
            m_infoCluster[rootW].parent = rootV;
            itName->second = rootV;
          }
        }
      }
    }
  }

  if (considerSmtDependencies && residual && !residual->empty()) {
    // Merge components connected via SMT assertions linking different theory
    // variables. Names without a corresponding CNF variable receive temporary
    // union-find nodes so they can bridge Boolean components.
    Var nextTmp = m_nbVar + m_clauses.size() + 1;

    if (traceDecompose) {
      Trace("decompose") << "c SMT assertions connecting variables:\n";
    }

    for (size_t assertionIdx = 0; assertionIdx < residual->size();
         assertionIdx++) {
      const auto &a = (*residual)[assertionIdx];
      const auto &vars = (*residualVarLists)[assertionIdx];
      if (vars.empty()) continue;

      // Print which SMT constraint connects which variables
      if (traceDecompose) {
        std::unordered_set<Var> cnfVarsInAssertion;
        for (const auto &name : vars) {
          if (residualVarNames->find(name) == residualVarNames->end()) continue;
          for (Var *it = m_activeVariables; it != lastVar; ++it) {
            Var v = *it;
            if (v >= (int)d4_smt_idxToTerm.size() ||
                m_currentValue[v] != l_Undef)
              continue;
            const auto &cnfVars = smtVarNamesCache[v];
            if (std::find(cnfVars.begin(), cnfVars.end(), name) !=
                cnfVars.end()) {
              cnfVarsInAssertion.insert(v);
            }
          }
        }

        std::ostringstream ss;
        ss << "c Assertion " << assertionIdx << " (" << a
           << ") connects:";

        ss << "\n  SMT vars: {";
        bool first = true;
        for (const auto &name : vars) {
          if (!first) ss << ", ";
          ss << name;
          first = false;
        }
        ss << "}";

        ss << "\n  CNF vars: {";
        first = true;
        for (const auto &v : cnfVarsInAssertion) {
          if (!first) ss << ", ";
          ss << v;
          first = false;
        }
        ss << "}";

        Trace("decompose") << ss.str() << "\n";
      }

      Var rep = -1;
      for (const auto &name : vars) {
        if (residualVarNames->find(name) == residualVarNames->end()) continue;
        auto itName = nameToRoot.find(name);
        if (itName == nameToRoot.end()) {
          Var idx = nextTmp++;
          if ((size_t)idx >= m_infoCluster.size())
            m_infoCluster.resize(idx + 1, {0, 0, -1});
          m_infoCluster[idx].parent = idx;
          m_infoCluster[idx].size = 1;
          itName = nameToRoot.emplace(name, idx).first;
        }

        Var rootW = findRoot(itName->second);

        if (rep == -1) {
          rep = rootW;
        } else if (rep != rootW) {
          if (traceDecompose) {
            auto itCompRep = componentTrace.find(rep);
            auto itCompW = componentTrace.find(rootW);
            if (itCompRep != componentTrace.end() && !itCompRep->second.empty() &&
                itCompW != componentTrace.end() && !itCompW->second.empty()) {
              Trace("decompose") << "Assertion " << assertionIdx
                                 << " with SMT variable '" << name
                                 << "' merges CNF components "
                                 << describeComponent(rep) << " and "
                                 << describeComponent(rootW) << std::endl;
            }
          }
          if (m_infoCluster[rep].size < m_infoCluster[rootW].size) {
            if (traceDecompose) mergeComponentData(rootW, rep);
            m_infoCluster[rootW].size += m_infoCluster[rep].size;
            m_infoCluster[rep].parent = rootW;
            rep = rootW;
          } else {
            if (traceDecompose) mergeComponentData(rep, rootW);
            m_infoCluster[rep].size += m_infoCluster[rootW].size;
            m_infoCluster[rootW].parent = rep;
          }
        }

        itName->second = rep;
      }
    }
  }

  // Print final components after merging
  if (traceDecompose) {
    // Group variables by their final component roots
    std::unordered_map<Var, std::vector<Var>> finalComponents;
    int finalComponentCount = 0;

    for (Var *it = m_activeVariables; it != lastVar; ++it) {
      Var v = *it;
      if (m_currentValue[v] != l_Undef) continue;

      // Find the root component
      Var rootV = findRoot(v);

      // Add variable to its component
      finalComponents[rootV].push_back(v);
    }

    // Count and print non-trivial components
    for (const auto& comp : finalComponents) {
      if (comp.second.size() > 1) finalComponentCount++;
    }

    Trace("decompose") << "c Final components after merging: "
                        << finalComponentCount << "\n";
    int compIdx = 0;
    for (const auto& comp : finalComponents) {
      if (comp.second.size() > 1) {
        std::ostringstream ss;
        ss << "c Final component " << compIdx++ << " vars: {";
        for (size_t i = 0; i < comp.second.size(); i++) {
          ss << comp.second[i];
          if (i < comp.second.size() - 1) ss << ", ";
        }
        ss << "}";
        Trace("decompose") << ss.str() << "\n";
      }
    }
  }

  if (traceDecompose) {
    int smtComponentCount = 0;
    if (d4_smt_solver) {
      // Group variables by component for printing
      std::unordered_map<Var, std::unordered_set<std::string>> componentToSmtVars;

      // First collect the components by their roots
      for (Var *it = m_activeVariables; it != lastVar; ++it) {
        Var v = *it;
        if (v >= (int)d4_smt_idxToTerm.size() || m_currentValue[v] != l_Undef) continue;
        const cvc5::Term &t = d4_smt_idxToTerm[v];
        if (t.isNull()) continue;

        // Find component root
        Var rootV = findRoot(v);

        // Add SMT variables to this component
        const auto &vars = smtVarNamesCache[v];
        componentToSmtVars[rootV].insert(vars.begin(), vars.end());
      }

      // Print components with their SMT variables
      if (traceDecompose) {
        Trace("decompose") << "c SMT variable components:\n";
        int compIdx = 0;
        for (const auto& comp : componentToSmtVars) {
          if (comp.second.empty()) continue;
          std::ostringstream ss;
          ss << "c Component " << compIdx++ << " SMT vars: {";
          bool first = true;
          for (const auto& name : comp.second) {
            if (!first) ss << ", ";
            ss << name;
            first = false;
          }
          ss << "}";
          Trace("decompose") << ss.str() << "\n";
        }
      }

      smtComponentCount = componentToSmtVars.size();
    }

    std::cout << "c [CONNECTED COMPONENT] CNF: " << cnfComponentCount
              << ", SMT: " << smtComponentCount << '\n';
  }

  // collect the component.
  std::vector<Var> rootSet;
  freeVar.resize(0);

  for (currentVar = m_activeVariables; currentVar != lastVar; currentVar++) {
    Var v = *currentVar;
    assert(m_currentValue[v] == l_Undef);

    if (getNbClause(v) == 0) {
      freeVar.push_back(v);
      continue;
    }
    if (m_infoCluster[v].parent == v && m_infoCluster[v].size == 1) {
      freeVar.push_back(v);
      continue;
    }
    assert(m_currentValue[v] == l_Undef);

    // get the root.
    unsigned rootV = findRoot(v);

    if (m_infoCluster[rootV].pos == -1) {
      m_infoCluster[rootV].pos = varCo.size();
      varCo.push_back(std::vector<Var>());
      rootSet.push_back(rootV);
    }

    varCo[m_infoCluster[rootV].pos].push_back(v);
  }

  // restore for the next run.
  resetUnMark();
  for (auto &v : rootSet) m_infoCluster[v].pos = -1;
  int nbComp = varCo.size();
  if (nbComp > 1) {
    Profile.addDecomposition();
  } else if (considerSmtDependencies && cnfComponentCount > 1) {
    Profile.addStoppedDecomposition();
    if (traceDecompose) {
      Trace("decompose")
          << "SMT dependencies merged " << cnfComponentCount
          << " candidate CNF components into a single block" << std::endl;
    }
  }
  return nbComp;
}  // computeConnectedComponent

/**
   Collect the set of literals connected to l and store the result in
   varComponent.

   @param[in] l, the considered literal
   @param[in] v, the label of the previously assigned component (0 if not
   assigned).
   @param[in] varComponent, the set of varaible connected to l.
   @param[in] nbComponent, the component label.
*/
void CnfManager::connectedToLit(Lit l, std::vector<int> &v,
                                std::vector<Var> &varComponent,
                                int nbComponent) {
  for (unsigned i = 0; i < 2; i++) {
    IteratorIdxClause listIndex =
        i ? getVecIdxClauseBin(l) : getVecIdxClauseNotBin(l);

    for (int *ptr = listIndex.start; ptr != listIndex.end; ptr++) {
      int idx = *ptr;

      // As in computeConnectedComponent, clause indices originating from the
      // occurrence list may become invalid.  Protect accesses to m_markView to
      // avoid out-of-bounds reads on corrupted indices.
      if (idx < 0 || static_cast<size_t>(idx) >= m_markView.size()) {
        continue;
      }

      if (m_markView[idx]) continue;
      m_markView[idx] = 1;
      m_mustUnMark.push_back(idx);

      // compute component
      for (auto &l : m_clauses[idx]) {
        if (m_currentValue[l.var()] != l_Undef || v[l.var()]) continue;

        varComponent.push_back(l.var());
        v[l.var()] = nbComponent;
      }
    }
  }
}  // connectedToLit

/**
   Look all the formula in order to compute the connected component
   of the formula (union find algorithm).

   @param[out] varCo, the different connected components found
   @param[in] setOfVar, the current set of variables
   @param[in] isProjected, a boolean vbector that spectify the targeted
   variables.
   @param[out] freeVar, the set of variables that are present in setOfVar but
   not in the problem anymore

   \return the number of component found
*/
int CnfManager::computeConnectedComponentTargeted(
    std::vector<std::vector<Var>> &varCo, std::vector<Var> &setOfVar,
    std::vector<bool> &isProjected, std::vector<Var> &freeVar) {
  freeVar.clear();
  std::vector<Var> projVar;
  projVar.reserve(setOfVar.size());
  for (Var v : setOfVar) {
    if (m_currentValue[v] != l_Undef) continue;
    if (isProjected[v])
      projVar.push_back(v);
    else
      freeVar.push_back(v);
  }

  std::vector<Var> extraFree;
  int nbComponent = computeConnectedComponent(varCo, projVar, extraFree);
  freeVar.insert(freeVar.end(), extraFree.begin(), extraFree.end());
  return nbComponent;
}  // computeConnectedComponentTargeted

/**
   Test if a given clause is actually satisfied under the current
   interpretation.

   @param[in] idx, the clause index.

   \return true if the clause is satisfied, false otherwise.
*/
bool CnfManager::isSatisfiedClause(unsigned idx) {
  assert(idx < m_clauses.size());
  return m_infoClauses[idx].isSat;
}  // isSatisfiedClause

/**
   Test if a given clause is actually satisfied under the current
   interpretation.

   @param[in] idx, the clause index.

   \return true if the clause is satisfied, false otherwise.
*/
bool CnfManager::isSatisfiedClause(std::vector<Lit> &c) {
  for (auto &l : c) {
    if (!litIsAssigned(l)) continue;
    if (l.sign() && m_currentValue[l.var()] == l_False) return true;
    if (!l.sign() && m_currentValue[l.var()] == l_True) return true;
  }

  return false;
}  // isSatisfiedClause

/**
   Test at the same time if a given clause is actually satisfied under
   the current interpretation and if its set of variables belong to the
   current component that is represented as a boolean map given in
   parameter.

   @param[in] idx, the clause index.
   @param[in] currentComponent, currentComponent[var] is true when var is in
   the current component, false otherwise.

   \return false if the clause is satisfied, true otherwise.
*/
bool CnfManager::isNotSatisfiedClauseAndInComponent(
    int idx, std::vector<bool> &inCurrentComponent) {
  if (m_infoClauses[idx].isSat) return false;
  assert(!litIsAssigned(m_clauses[idx][0]));
  return inCurrentComponent[m_clauses[idx][0].var()];
}  // isSatisfiedClause

void CnfManager::showFormula(std::ostream &out) {
  out << "p cnf " << getNbVariable() << " " << getNbClause() << "\n";
  for (auto &cl : m_clauses) {
    showListLit(out, cl);
    out << "0\n";
  }
}  // showFormula

/**
 * @brief CnfManager::showCurrentFormula implementation.
 */
void CnfManager::showCurrentFormula(std::ostream &out) {
  out << "p cnf " << getNbVariable() << " " << getNbClause() << "\n";
  for (unsigned i = 0; i < m_clauses.size(); i++) {
    if (m_infoClauses[i].isSat) continue;
    out << "[" << i << "] ";
    for (auto &l : m_clauses[i])
      if (!litIsAssigned(l)) out << l << " ";
    out << "0\n";
  }
}  // showFormula

/**
 * @brief CnfManager::showCurrentFormula implementation.
 */
void CnfManager::showCurrentFormula(std::ostream &out,
                                    std::vector<bool> &isInComponent) {
  out << "p cnf " << getNbVariable() << " " << getNbClause() << "\n";
  for (unsigned i = 0; i < m_clauses.size(); i++) {
    if (!isNotSatisfiedClauseAndInComponent(i, isInComponent)) continue;
    if (m_infoClauses[i].isSat) continue;
    for (auto &l : m_clauses[i])
      if (!litIsAssigned(l)) out << l << " ";
    out << "0\n";
  }
}  // showFormula

/**
 * @brief CnfManager::debugFunction implementation.
 */
void CnfManager::debugFunction() {
  for (unsigned i = 0; i < m_clauses.size(); i++) {
    if (m_infoClauses[i].isSat) continue;
    for (auto &l : m_clauses[i])
      if (!litIsAssigned(l)) {
        unsigned nbOcc = 0;
        for (unsigned j = 0; j < m_occurrence[l.intern()].nbBin; j++)
          if (i == m_occurrence[l.intern()].bin[j]) nbOcc++;
        for (unsigned j = 0; j < m_occurrence[l.intern()].nbNotBin; j++)
          if (i == m_occurrence[l.intern()].notBin[j]) nbOcc++;

        if (nbOcc != 1) {
          std::cout << "literal " << l << " nbOcc = " << nbOcc
                    << " and the clause is " << i << '\n';
        }
        assert(nbOcc == 1);
      }
  }

  for (unsigned i = 0; i < m_occurrence.size(); i++) {
    if (m_currentValue[i >> 1] != l_Undef) continue;
    if (!m_occurrence[i].nbBin && !m_occurrence[i].nbNotBin) continue;
    for (unsigned j = 0; j < m_occurrence[i].nbBin; j++)
      assert(!m_infoClauses[m_occurrence[i].bin[j]].isSat);
    for (unsigned j = 0; j < m_occurrence[i].nbNotBin; j++) {
      if (m_infoClauses[m_occurrence[i].notBin[j]].isSat)
        std::cout << "=> " << (i >> 1) << " " << m_occurrence[i].notBin[j]
                  << '\n';
      assert(!m_infoClauses[m_occurrence[i].notBin[j]].isSat);
    }
  }

}  // debugFunction

}  // namespace d4
