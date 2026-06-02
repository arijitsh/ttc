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

#include "WrapperMinisat.hpp"

#include <cassert>
#include <iostream>
#include <typeinfo>
#include "../../../../ddnnf/logger.hpp"
#include "../../../../ddnnf/profiler.hpp"
#include "../../../../features.hpp"

#include "minisat/Solver.hpp"
#include "minisat/SolverTypes.hpp"
#include "minisat/mtl/Vec.hpp"
#include "src/problem/CnfMatrix.hpp"
#include "src/problem/ProblemManager.hpp"
#include "src/utils/ErrorCode.hpp"

namespace d4 {
using minisat::toInt;

/**
   This function initializes the SAT solver with a given problem.  Warning: we
   suppose that p is a CNF, otherwise a bad_cast exception is threw.

   @param[in] p, the problem we want to link with the SAT solver.
 */
void WrapperMinisat::initSolver(ProblemManager &p) {
  try {
    CnfMatrix &pcnf = dynamic_cast<CnfMatrix &>(p);

    // say to the solver we have pcnf.getNbVar() variables.
    while ((unsigned)m_solver.nVars() <= p.getNbVar()) m_solver.newVar();
    m_model.resize(p.getNbVar() + 1, l_Undef);

    // load the clauses
    std::vector<std::vector<Lit>> &clauses = pcnf.getClauses();
    for (auto &cl : clauses) {
      minisat::vec<minisat::Lit> lits;
      for (auto &l : cl) lits.push(minisat::mkLit(l.var(), l.sign()));
      m_solver.addClause(lits);
    }
  } catch (std::bad_cast &bc) {
    std::cerr << "c bad_cast caught: " << bc.what() << '\n';
    std::cerr << "c A CNF formula was expeted\n";
    exit(ERROR_BAD_CAST);
  }

  m_activeModel = false;
  m_needModel = false;
  setNeedModel(m_needModel);
  m_isInAssumption.resize(p.getNbVar() + 1, 0);
  m_syncedTrailLits.clear();
  m_smtTrail = 0;
  m_smtContextLevel = 0;
  m_pendingTheoryCheck = false;
}  // initSolver

/**
   Call the SAT solver and return its result.

   @param[in] setOfvar, the variables we focus the solving process.

   \return true if the problem is SAT, false otherwise.
 */
bool WrapperMinisat::solve(std::vector<Var> &setOfVar) {
  if (Trace.isEnabled("d4")) {
    Trace("d4") << "c [WrapperMinisat] solve with connected component of size "
                 << setOfVar.size() << "\n";
  }
  // Record variables of the current connected component when
  // component-based decomposition is enabled.
  if (!d4_no_decompose)
    m_currentComponent.assign(setOfVar.begin(), setOfVar.end());
  if (m_activeModel && m_needModel) {
    // Ensure the SMT solver is aware of the current Boolean assignment and
    // that the assignment satisfies theory constraints when all projection
    // variables have been instantiated.
    syncTrail();
    bool smtSat = true;
    if (d4_smt_solver && allProjAssigned()) {
      Profile.addCheckSatCall();
      cvc5::Result r = d4_smt_solver->checkSat();
      if (r.isUnsat()) Profile.addUnsatCheckSatCall();
      smtSat = r.isSat();
      m_pendingTheoryCheck = false;
      if (!smtSat) m_activeModel = false;
    }
    if (Trace.isEnabled("d4")) {
      auto stream = Trace("d4");
      stream << "trail:";
      for (int i = 0; i < m_solver.trail.size(); ++i) {
        minisat::Lit l = m_solver.trail[i];
        stream << ' ' << (minisat::sign(l) ? '-' : '+') << minisat::var(l);
      }
      stream << " -> CNF SAT";
      if (d4_smt_solver && allProjAssigned())
        stream << ", SMT " << (smtSat ? "SAT" : "UNSAT");
      stream << std::endl;
    }
    m_lastCnfSat = true;
    m_lastSmtSat = smtSat;
    return m_activeModel;
  }

  m_setOfVar_m.setSize(0);
  for (auto &v : setOfVar) m_setOfVar_m.push(v);
  double rebuildStart = Log.elapsed();
  m_solver.rebuildWithConnectedComponent(m_setOfVar_m);
  double rebuildEnd = Log.elapsed();

  Profile.addSatCall();
  bool cnfSat = m_solver.solveWithAssumptions();
  if (!cnfSat) Profile.addUnsatSatCall();
  syncTrail();
  bool smtSat = true;
  if (cnfSat && d4_smt_solver && allProjAssigned()) {
    Profile.addCheckSatCall();
    cvc5::Result r = d4_smt_solver->checkSat();
    if (r.isUnsat()) Profile.addUnsatCheckSatCall();
    smtSat = r.isSat();
    m_pendingTheoryCheck = false;
  }
  if (Trace.isEnabled("d4")) {
    auto stream = Trace("d4");
    stream << "trail:";
    for (int i = 0; i < m_solver.trail.size(); ++i) {
      minisat::Lit l = m_solver.trail[i];
      stream << ' ' << (minisat::sign(l) ? '-' : '+') << minisat::var(l);
    }
    stream << " -> CNF " << (cnfSat ? "SAT" : "UNSAT");
    if (d4_smt_solver && allProjAssigned())
      stream << ", SMT " << (smtSat ? "SAT" : "UNSAT");
    stream << std::endl;
  }
  m_lastCnfSat = cnfSat;
  m_lastSmtSat = smtSat;
  m_activeModel = cnfSat && smtSat;
  print_progress_d4(0.0);
  return m_activeModel;
}  // solve

/**
   Call the SAT solver and return its result.

   \return true if the problem is SAT, false otherwise.
 */
bool WrapperMinisat::solve() {
  std::cout << "c [WrapperMinisat] solve with all variables\n";
  // The entire formula is considered a single component here. Avoid the
  // bookkeeping overhead when decomposition is disabled.
  if (!d4_no_decompose) {
    m_currentComponent.clear();
    for (int v = 0; v < m_solver.nVars(); ++v) m_currentComponent.push_back(v);
  }
  double rebuildStart = Log.elapsed();
  m_solver.rebuildWithAllVar();
  double rebuildEnd = Log.elapsed();
  Profile.addSatCall();
  bool cnfSat = m_solver.solveWithAssumptions();
  if (!cnfSat) Profile.addUnsatSatCall();
  syncTrail();
  bool smtSat = true;
  if (cnfSat && d4_smt_solver && allProjAssigned()) {
    Profile.addCheckSatCall();
    cvc5::Result r = d4_smt_solver->checkSat();
    if (r.isUnsat()) Profile.addUnsatCheckSatCall();
    smtSat = r.isSat();
    m_pendingTheoryCheck = false;
  }
  if (Trace.isEnabled("d4")) {
    auto stream = Trace("d4");
    stream << "trail:";
    for (int i = 0; i < m_solver.trail.size(); ++i) {
      minisat::Lit l = m_solver.trail[i];
      stream << ' ' << (minisat::sign(l) ? '-' : '+') << minisat::var(l);
    }
    stream << " -> CNF " << (cnfSat ? "SAT" : "UNSAT");
    if (d4_smt_solver && allProjAssigned())
      stream << ", SMT " << (smtSat ? "SAT" : "UNSAT");
    stream << std::endl;
  }
  m_lastCnfSat = cnfSat;
  m_lastSmtSat = smtSat;
  print_progress_d4(0.0);
  return cnfSat && smtSat;
}  // solve

/**
 * @brief Enforce the unit propagation of all the assumption literals.
 *
 * @return true if we did not reach a conflict, false otherwise.
 */
bool WrapperMinisat::propagateAssumption() {
  bool ok = m_solver.propagateAssumption();
  syncTrail();
  return ok;
}  // propagateAssumption

void WrapperMinisat::addClause(const std::vector<Lit> &clause) {
  // Ensure clause injection occurs at the top level. Minisat expects
  // decisionLevel() == 0 when adding a clause, so we explicitly
  // backtrack before handing the clause to the solver.
  m_solver.cancelUntil(0);
  minisat::vec<minisat::Lit> lits;
  for (const auto &l : clause) lits.push(minisat::mkLit(l.var(), l.sign()));
  m_solver.addClause(lits);
}

/**
 * @brief Strong assumption in the sense we push the literal on the stack.
 *
 * @param l the literal we want to push.
 */
void WrapperMinisat::uncheckedEnqueue(Lit l) {
  m_solver.uncheckedEnqueue(minisat::mkLit(l.var(), l.sign()));
}  // uncheckedEnqueue

/**
   An accessor on the activity of a variable.

   @param[in] v, the variable we want the activity.

   \return the activity of v.
 */
double WrapperMinisat::getActivity(Var v) {
  return m_solver.activity[v];
}  // getActivity

/**
 * @brief Set the reverse polarity flag to the solver.
 *
 * @param value is the value we want to assign.
 */
void WrapperMinisat::setReversePolarity(bool value) {
  m_solver.reversePolarity = value;
}  // setReversePolarity

/**
 * @brief Return the number of times the variable v occurs in a conflict.
 * @param[in] v, the variable we want the activity.
 * \return the number of times v occurs in a conflict.
 */
double WrapperMinisat::getCountConflict(Var v) {
  assert(v >= 0);
  assert(v < m_solver.scoreActivity.size());
  return m_solver.scoreActivity[v];
}  // getCountConflict

/**
 * @brief decayCountConflict implementation.
 *
 */
void WrapperMinisat::decayCountConflict() {
  for (unsigned i = 0; i < m_solver.scoreActivity.size(); i++)
    m_solver.scoreActivity[i] = m_solver.scoreActivity[i] / 2;
}  // decayCountConflict

/**
 * @brief Set the count conflict in the solver.
 * @param[in] v is the variable we want to set the counter of  conflicts.
 * @param[in] count is the count we want to assign.
 */
void WrapperMinisat::setCountConflict(Var v, double count) {
  m_solver.scoreActivity[v] = count;
}  // setCountConflict

/**
   Print out the trail on the standard output.
 */
void WrapperMinisat::showTrail() { m_solver.showTrail(); }  // showTrail

/**
   An accessor on the polarity of a variable.

   @param[in] v, the variable we want the polarity.
 */
bool WrapperMinisat::getPolarity(Var v) {
  return m_solver.polarity[v];
}  // getPolarity

/**
   Collect the unit literal from the affectation of the literal l to the
   formula.

   @param[in] l, the literal we want to branch on.
   @param[out] units, the unit literals

   \return true if assigning l and propagating it does not give a conflict,
   false otherwise.
 */
bool WrapperMinisat::decideAndComputeUnit(Lit l, std::vector<Lit> &units) {
  if (!m_solver.okay()) return false;
  minisat::Lit ml = minisat::mkLit(l.var(), l.sign());
  if (varIsAssigned(l.var())) {
    if (m_solver.litAssigned(l.var()) != ml) return false;
    units.push_back(l);
    return true;
  }

  int posTrail = (m_solver.trail).size();
  m_solver.newDecisionLevel();
  m_solver.uncheckedEnqueue(ml);
  minisat::CRef confl = m_solver.propagate();
  syncTrail();

  if (confl != minisat::CRef_Undef)  // unit literal
  {
    int bt;
    minisat::vec<minisat::Lit> learnt_clause;
    m_solver.analyzeLastUIP(confl, learnt_clause, bt);
    m_solver.cancelUntil(m_solver.decisionLevel() - 1);
    syncTrail();
    assert(learnt_clause[0] == minisat::mkLit(l.var(), !l.sign()));
    m_solver.insertClauseAndPropagate(learnt_clause);
    return false;
  }

  for (int j = posTrail; j < m_solver.trail.size(); j++)
    units.push_back(
        Lit::makeLit(var(m_solver.trail[j]), sign(m_solver.trail[j])));
  m_solver.cancelUntil(m_solver.decisionLevel() - 1);
  syncTrail();
  return true;
}  // decideAndComputeUnit

/**
 * @brief WrapperMinisat::literalProbing implementation.
 */
bool WrapperMinisat::failedLiteralProbing(Lit l) {
  if (!m_solver.okay()) return true;
  minisat::Lit ml = minisat::mkLit(l.var(), (~l).sign());
  if (varIsAssigned(l.var())) {
    if (m_solver.litAssigned(l.var()) == ml) return true;
    return false;
  }

  m_solver.newDecisionLevel();
  m_solver.uncheckedEnqueue(ml);
  minisat::CRef confl = m_solver.propagate();
  syncTrail();
  m_solver.cancelUntil(m_solver.decisionLevel() - 1);
  syncTrail();

  if (confl != minisat::CRef_Undef) return true;  // unit literal
  return false;
}  // failedLiteralProbing

/**
   Fill the vector units with the literal l that are units such that l.var() is
   in component.

   @param[in] component, the set of variables we search for.
   @param[out] units, the place where we store the literals found.
 */
void WrapperMinisat::whichAreUnits(std::vector<Var> &component,
                                   std::vector<Lit> &units) {
  auto add_unique = [&](Lit l) {
    for (const auto &u : units) {
      if (u.var() == l.var()) return;  // already present
    }
    units.push_back(l);
  };

  // Add all variables assigned by the SAT solver within the current component.
  for (auto &v : component) {
    if (!m_solver.isAssigned(v)) continue;
    minisat::Lit l = m_solver.litAssigned(v);
    add_unique(Lit::makeLit(var(l), sign(l)));
  }

  // Also include explicit assumptions so that the CNF manager sees them
  // as fixed at this recursion level, even if the underlying SAT solver
  // keeps them at non-zero decision levels.
  if (!m_assumption.empty()) {
    for (const auto &l : m_assumption) {
      Var v = l.var();
      // Include only if the variable is part of this component.
      for (auto &cv : component) {
        if (cv == v) {
          add_unique(l);
          break;
        }
      }
    }
  }
}  // whichAreUnits

/**
 * @brief Get the list of unit literals that are in the trail (we suppose that
 * the decision level is zero).
 *
 * @param[out] units is the list of unit literals.
 */
void WrapperMinisat::getUnits(std::vector<Lit> &units) {
  for (int i = 0; i < m_solver.trail.size(); i++) {
    minisat::Lit l = m_solver.trail[i];
    units.push_back(Lit::makeLit(var(l), sign(l)));
  }
}  // getUnits

/**
   Check out if the given variable is assigned or not by the solver.

   @param[in] v, the variable we search for.

   \return true if the variable is assigned, false otherwise.
 */
bool WrapperMinisat::varIsAssigned(Var v) {
  return m_solver.isAssigned(v);
}  // varIsAssigned

/**
   Restart the solver.
 */
void WrapperMinisat::restart() {
  m_solver.cancelUntil(0);
  syncTrail();
}  // restart

/**
   Transfer to the solver the fact we have a set of assumption variables we want
   to consider.

   @param[in] assums, the set of assumptions
 */
void WrapperMinisat::setAssumption(std::vector<Lit> &assums) {
  popAssumption(m_assumption.size());
  minisat::vec<minisat::Lit> &assumptions = m_solver.assumptions;
  assumptions.clear();
  m_assumption.clear();
  for (auto &l : assums) pushAssumption(l);
  syncTrail();
}  // setAssumption

/**
   \return the current assumption.

   @param[in] assums, the set of assumptions
 */
std::vector<Lit> &WrapperMinisat::getAssumption() {
  return m_assumption;
}  // getAssumption

/**
   Print out the assumption.

   @param[in] out, the stream where is print the assumption.
 */
void WrapperMinisat::displayAssumption(std::ostream &out) {
  minisat::vec<minisat::Lit> &assumptions = m_solver.assumptions;
  for (int i = 0; i < assumptions.size(); i++) {
    minisat::Lit l = assumptions[i];
    std::cout << (minisat::sign(l) ? "-" : "") << minisat::var(l) << " ";
  }
  std::cout << "\n";
}  // displayAssumption

/**
   Ask for the model.

   @param[in] b, a boolean value to ask the solver to get the model.
 */
void WrapperMinisat::setNeedModel(bool b) {
  m_needModel = b;
  m_solver.setNeedModel(b);
}  // setNeedModel

/**
 * @brief Return the model computed by the solver.
 *
 * @return the model's value (lbool).
 */
std::vector<lbool> &WrapperMinisat::getModel() {
  for (int i = 0; i < m_solver.model.size(); i++) {
    if (minisat::toInt(m_solver.model[i]) == 0)
      m_model[i] = l_True;
    else if (minisat::toInt(m_solver.model[i]) == 1)
      m_model[i] = l_False;
    else
      m_model[i] = l_Undef;
  }

  return m_model;
}  // getModel

/**
 * @brief Get the value given by the last computed model.
 *
 * @param v is the variable we want to get the assignment.
 * @return the last value of v.
 */
lbool WrapperMinisat::getModelVar(Var v) {
  return minisat::toInt(m_solver.model[v]);
}  // getModelVar

/**
   Push a new assumption.

   @param[in] l, the literal we want to push.
 */
void WrapperMinisat::pushAssumption(Lit l) {
  minisat::Lit ml = minisat::mkLit(l.var(), l.sign());

  if (m_isInAssumption[l.var()]) {
    if ((m_isInAssumption[l.var()] - 1) == l.sign()) return;
    for (unsigned i = 0; i < m_assumption.size(); ++i) {
      if (m_assumption[i].var() == l.var()) {
        (m_solver.assumptions)[i] = ml;
        m_assumption[i] = l;
        // A changed assumption invalidates any previously stored model.
        m_activeModel = false;
        m_isInAssumption[l.var()] = 1 + l.sign();
        if (false && m_activeModel && m_needModel) {
          m_activeModel = m_solver.litTrueInLastModel(ml);
          if (m_activeModel) {
            assert(m_solver.decisionLevel() == m_solver.assumptions.size() - 1);
            m_solver.newDecisionLevel();
            assert(!m_solver.isAssigned(var(ml)));
            m_solver.uncheckedEnqueue(ml);
            [[maybe_unused]] minisat::CRef cref = m_solver.propagate();
            syncTrail();
            assert(cref == minisat::CRef_Undef);
            m_solver.cancelUntil(m_solver.decisionLevel() - 1);
            syncTrail();
          }
        }
        return;
      }
    }
  }

  // Introducing a new assumption requires a fresh SAT check.
  m_activeModel = false;
  (m_solver.assumptions).push(ml);
  m_assumption.push_back(l);
  m_isInAssumption[l.var()] = 1 + l.sign();

  if (false && m_activeModel && m_needModel) {
    m_activeModel = m_solver.litTrueInLastModel(ml);
    if (m_activeModel) {
      assert(m_solver.decisionLevel() == m_solver.assumptions.size() - 1);
      m_solver.newDecisionLevel();
      assert(!m_solver.isAssigned(var(ml)));
      m_solver.uncheckedEnqueue(ml);
      [[maybe_unused]] minisat::CRef cref = m_solver.propagate();
      syncTrail();
      assert(cref == minisat::CRef_Undef);
      m_solver.cancelUntil(m_solver.decisionLevel() - 1);
      syncTrail();
    }
  }
}  // pushAssumption

/**
   Remove the last assumption and cancelUntil.

   @param[in] count, the number of element we pop.
 */
void WrapperMinisat::popAssumption(unsigned count) {
  for (unsigned i = m_assumption.size() - count; i < m_assumption.size(); i++) {
    assert(m_isInAssumption[m_assumption[i].var()]);
    m_isInAssumption[m_assumption[i].var()] = 0;
  }

  m_assumption.resize(m_assumption.size() - count);
  (m_solver.assumptions).shrink_(count);
  (m_solver.cancelUntil)((m_solver.assumptions).size());
  syncTrail();
}  // popAssumption

inline unsigned WrapperMinisat::getNbConflict() { return m_solver.conflicts; }
bool WrapperMinisat::isUnsat() { return !m_solver.okay(); }

/**
 * @brief Compute the core.
 *
 */
void WrapperMinisat::getCore() {
  for (unsigned i = 0; i < m_solver.conflict.size(); i++) {
    minisat::Lit l = m_solver.conflict[i];
    std::cout << (minisat::sign(l) ? "-" : "") << minisat::var(l) << "("
              << m_solver.level(var(l)) << ") ";
  }
  std::cout << "  ---> " << m_solver.decisionLevel() << "\n";
}  // getCore

/**
 * @brief
 *
 * @param l
 */
void WrapperMinisat::getLastIUP(Lit dl) {
  minisat::Lit ml = minisat::mkLit(dl.var(), dl.sign());
  if (m_solver.reason(minisat::var(ml)) == minisat::CRef_Undef) {
    std::cout << "decision\n";
    return;
  }

  minisat::vec<minisat::Lit> conf;
  m_solver.analyzeFinal(ml, conf);

  for (unsigned i = 0; i < conf.size(); i++) {
    minisat::Lit l = conf[i];
    std::cout << (minisat::sign(l) ? "-" : "") << minisat::var(l) << "("
              << m_solver.level(var(l)) << ") ";
  }
  std::cout << "  ---> " << m_solver.decisionLevel() << "\n";
}  // getLastIUP

void WrapperMinisat::syncTrail() {
  if (!d4_smt_solver) {
    m_syncedTrailLits.clear();
    m_smtTrail = 0;
    m_smtContextLevel = 0;
    m_pendingTheoryCheck = false;
    return;
  }

  auto& tm = ttc::getTermBuilder(*d4_smt_solver);

  const int trailSize = m_solver.trail.size();
  const int syncedSize = m_syncedTrailLits.size();
  int prefix = 0;
  while (prefix < syncedSize && prefix < trailSize) {
    minisat::Lit recorded = m_syncedTrailLits[prefix];
    minisat::Lit current = m_solver.trail[prefix];
    if (minisat::var(recorded) == minisat::var(current) &&
        minisat::sign(recorded) == minisat::sign(current)) {
      ++prefix;
    } else {
      break;
    }
  }

  if (prefix < syncedSize) {
    int pops = 0;
    for (int idx = syncedSize - 1; idx >= prefix; --idx) {
      minisat::Lit lit = m_syncedTrailLits.back();
      int v = minisat::var(lit);
      if (v < (int)d4_smt_idxToTerm.size()) {
        const cvc5::Term &term = d4_smt_idxToTerm[v];
        if (!term.isNull()) ++pops;
      }
      m_syncedTrailLits.pop_back();
    }
    if (pops > 0) {
      d4_smt_solver->pop(pops);
      m_smtContextLevel -= pops;
      if (m_smtContextLevel < 0) m_smtContextLevel = 0;
      m_pendingTheoryCheck = false;
    }
  }

  for (int i = prefix; i < trailSize; ++i) {
    minisat::Lit l = m_solver.trail[i];
    m_syncedTrailLits.push_back(l);
    int v = minisat::var(l);
    if (v >= (int)d4_smt_idxToTerm.size()) continue;
    cvc5::Term t = d4_smt_idxToTerm[v];
    if (t.isNull()) continue;
    if (m_solver.reason(v) == minisat::CRef_Undef && m_solver.level(v) > 0) {
      // Decision counting now handled alongside SMT cache lookups.
    }
    if (Trace.isEnabled("d4") && m_solver.reason(v) == minisat::CRef_Undef &&
        m_solver.level(v) > 0) {
      auto decideStream = Trace("d4");
      decideStream << "decide " << (minisat::sign(l) ? '-' : '+') << v
                   << " lvl " << m_solver.level(v);
      if (!t.isNull()) decideStream << " (" << t << ")";
      decideStream << std::endl;

      auto cnfStream = Trace("d4");
      cnfStream << "CNF state:";
      for (int j = 0; j < i; ++j) {
        minisat::Lit lj = m_solver.trail[j];
        cnfStream << ' ' << (minisat::sign(lj) ? '-' : '+') << minisat::var(lj);
      }
      cnfStream << std::endl;

      auto assertions = d4_smt_solver->getAssertions();
      auto smtStream = Trace("d4");
      smtStream << "SMT state:";
      for (const auto &a : assertions) smtStream << ' ' << a;
      smtStream << std::endl;
    }
    d4_smt_solver->push();
    if (minisat::sign(l))
      d4_smt_solver->assertFormula(tm.mkTerm(cvc5::Kind::NOT, {t}));
    else
      d4_smt_solver->assertFormula(t);
    ++m_smtContextLevel;
    m_pendingTheoryCheck = true;
  }

  m_smtTrail = m_syncedTrailLits.size();
  print_progress_d4(0.0);

  if (Trace.isEnabled("trail")) {
    auto stream = Trace("trail");
    bool first = true;
    for (const auto &lit : m_syncedTrailLits) {
      int v = minisat::var(lit);
      if (v >= (int)d4_smt_idxToTerm.size()) continue;
      const cvc5::Term &term = d4_smt_idxToTerm[v];
      if (term.isNull()) continue;
      if (!first) stream << '\n';
      stream << ' ' << (minisat::sign(lit) ? '-' : '+') << v << '@'
             << m_solver.level(v) << '['
             << (m_solver.reason(v) == minisat::CRef_Undef ? 'd' : 'p')
             << "] -> " << term;
      first = false;
    }
    stream << " (depth " << m_smtContextLevel << ")" << std::endl;
  }

  if (Trace.isEnabled("trail-bool")) {
    auto stream = Trace("trail-bool");
    stream << "c Boolean trail:";
    const int trailSize = m_solver.trail.size();
    for (int i = 0; i < trailSize; ++i) {
      minisat::Lit lit = m_solver.trail[i];
      int v = minisat::var(lit);
      stream << ' ' << (minisat::sign(lit) ? '-' : '+') << v << '@'
             << m_solver.level(v) << '['
             << (m_solver.reason(v) == minisat::CRef_Undef ? 'd' : 'p')
             << ']';
      if (v < (int)d4_smt_idxToTerm.size()) {
        const cvc5::Term &term = d4_smt_idxToTerm[v];
        if (!term.isNull()) stream << " -> " << term;
      }
    }
    stream << std::endl;
  }
}

bool WrapperMinisat::allProjAssigned() {
  if (!d4_smt_solver) return false;

  int expectedTheory = 0;
  const int trailSize = m_solver.trail.size();
  for (int i = 0; i < trailSize; ++i) {
    minisat::Lit l = m_solver.trail[i];
    int v = minisat::var(l);
    if (v >= (int)d4_smt_idxToTerm.size()) continue;
    const cvc5::Term &term = d4_smt_idxToTerm[v];
    if (term.isNull()) continue;
    ++expectedTheory;
  }

  if (expectedTheory == 0) return false;
  if (expectedTheory != m_smtContextLevel) return false;
  return m_pendingTheoryCheck;
}

}  // namespace d4
